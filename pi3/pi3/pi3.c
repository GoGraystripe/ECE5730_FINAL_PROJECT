/* 2D dilute hard-disk gas, Pico 2 / RP2350.
 * Lab 2 VGA wiring: sync GPIO16/17, green18/19, blue20, red21.
 * Units: pixel, second, unit mass, k_B = 1. In two dimensions
 * K/N = k_B T and P*A = N*k_B*T (= nRT); P is force/length.
 * Finite-radius disks have excluded-volume corrections: the ideal law is
 * approached at low packing fraction, not enforced by rescaling velocities.
 * No thermostat, gravity, or forces between collisions.
 */
#include <math.h>
#include <stdint.h>
#include <string.h>
#ifndef GAS_HOST_TEST
#include "pico/stdlib.h"
#include "pico/multicore.h"
#include "pico/mutex.h"
#include "pico/rand.h"
#include "hardware/clocks.h"
#include "VGA/vga16_graphics_v3.h"
#endif

#ifndef GAS_PARTICLES
#define GAS_PARTICLES 80
#endif
#define GAS_RADIUS 2
// Reserve the left column for the project title and subtitle.
#define GAS_LEFT 224
#define GAS_RIGHT 615
#define GAS_TOP 24
#define GAS_BOTTOM 455
#define GAS_INITIAL_TEMPERATURE 1800.0f /* rms speed = sqrt(2*T) = 60 px/s */
#define GAS_STEP_US 4167u             /* fixed ~240 Hz, independent of VGA */
#define GAS_MASS 1.0f
#define GAS_BOLTZMANN 1.0f
#define GAS_PRESSURE_SECONDS 5.0f
#define GAS_DT (GAS_STEP_US * 1.e-6f)
#define XMIN (GAS_LEFT + GAS_RADIUS + 1.0f)
#define XMAX (GAS_RIGHT - GAS_RADIUS - 1.0f)
#define YMIN (GAS_TOP + GAS_RADIUS + 1.0f)
#define YMAX (GAS_BOTTOM - GAS_RADIUS - 1.0f)
#define WIDTH (GAS_RIGHT - GAS_LEFT - 2.0f)
#define HEIGHT (GAS_BOTTOM - GAS_TOP - 2.0f)
_Static_assert(GAS_PARTICLES > 1, "At least two particles required");
_Static_assert(GAS_RADIUS > 0 && GAS_LEFT >= 0 && GAS_RIGHT < 640 &&
               GAS_TOP >= 0 && GAS_BOTTOM < 480, "Invalid VGA geometry");
_Static_assert(XMAX > XMIN && YMAX > YMIN, "Container too small");

typedef struct { float x, y, vx, vy; } Particle;
static Particle particles[GAS_PARTICLES];
typedef struct {
    float kinetic_energy, temperature, pressure, ideal_pressure, compressibility;
    float pressure_x, pressure_y;
    uint32_t steps, pair_collisions, wall_collisions, late_steps, max_render_us;
} GasDiagnostics;
/* Inspect in debugger; only core 0 writes thermodynamic quantities.
 * Pressure is a non-overlapping simulated-time window, initially zero.
 * ideal_pressure and Z use the same window-averaged temperature as pressure. */
volatile GasDiagnostics gas_diagnostics;
static float impulse_x, impulse_y, pressure_time, energy_integral;
static uint32_t rng_state = 0x57301234u;
static uint32_t random_word(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}
static float uniform(void) { return ((random_word() >> 8) + 0.5f) / 16777216.0f; }
static float energy(void) {
    float k = 0;
    for (int i = 0; i < GAS_PARTICLES; ++i)
        k += .5f * GAS_MASS * (particles[i].vx * particles[i].vx +
                                     particles[i].vy * particles[i].vy);
    return k;
}
static int initialize_gas(void) {
    float mean_x = 0, mean_y = 0;
    for (int i = 0; i < GAS_PARTICLES; ++i) {
        int placed = 0;
        for (int attempt = 0; attempt < 20000 && !placed; ++attempt) {
            particles[i].x = XMIN + (XMAX-XMIN)*uniform();
            particles[i].y = YMIN + (YMAX-YMIN)*uniform();
            placed = 1;
            for (int j = 0; j < i; ++j) {
                float dx = particles[i].x-particles[j].x;
                float dy = particles[i].y-particles[j].y;
                if (dx*dx + dy*dy < (2*GAS_RADIUS+0.05f)*(2*GAS_RADIUS+0.05f)) {
                    placed = 0; break;
                }
            }
        }
        if (!placed) return 0; /* invalid/dense configuration: never spawn overlaps */
        float speed = sqrtf(-2.0f * logf(uniform()));
        float angle = 6.28318530718f * uniform();
        particles[i].vx = speed*cosf(angle);
        particles[i].vy = speed*sinf(angle);
        mean_x += particles[i].vx; mean_y += particles[i].vy;
    }
    for (int i = 0; i < GAS_PARTICLES; ++i) {
        particles[i].vx -= mean_x/GAS_PARTICLES;
        particles[i].vy -= mean_y/GAS_PARTICLES;
    }
    float scale = sqrtf(GAS_PARTICLES*GAS_BOLTZMANN*GAS_INITIAL_TEMPERATURE/energy());
    for (int i = 0; i < GAS_PARTICLES; ++i) {
        particles[i].vx *= scale; particles[i].vy *= scale;
    }
    return 1;
}

/* Earliest event within the remaining part of a fixed tick. Broad rejection
 * avoids sqrt for distant pairs. Swept tests prevent tunneling at any speed.
 * Equal-time events are resolved successively at zero elapsed time.
 * We never discard remaining time or cap velocity to meet a frame deadline. */
static void physics_step(void) {
    float remaining = GAS_DT;
    while (remaining > 0.0f) {
        float soonest = remaining;
        int hit_a = -1, hit_b = -1, wall_axis = -1;
        for (int i = 0; i < GAS_PARTICLES; ++i) {
            Particle *p = &particles[i];
            float t;
            if (p->vx != 0) {
                t = ((p->vx > 0 ? XMAX : XMIN)-p->x)/p->vx;
                if (t <= soonest) { soonest = fmaxf(t,0); hit_a=i; hit_b=-1; wall_axis=0; }
            }
            if (p->vy != 0) {
                t = ((p->vy > 0 ? YMAX : YMIN)-p->y)/p->vy;
                if (t <= soonest) { soonest = fmaxf(t,0); hit_a=i; hit_b=-1; wall_axis=1; }
            }
            for (int j = i+1; j < GAS_PARTICLES; ++j) {
                Particle *q = &particles[j];
                float dx=q->x-p->x, dy=q->y-p->y;
                float ux=q->vx-p->vx, uy=q->vy-p->vy;
                const float diameter = 2.0f*GAS_RADIUS;
                if (fabsf(dx) > diameter+fabsf(ux)*soonest ||
                    fabsf(dy) > diameter+fabsf(uy)*soonest) continue;
                float b=dx*ux+dy*uy;
                if (b >= 0) continue; /* separating or tangent */
                float a=ux*ux+uy*uy;
                float c=dx*dx+dy*dy-diameter*diameter;
                float disc=b*b-a*c;
                if (disc <= 0 || a == 0) continue;
                /* Stable smaller root; contact roundoff resolves immediately. */
                t = c > 0 ? c/(-b+sqrtf(disc)) : 0;
                if (t <= soonest) { soonest=t; hit_a=i; hit_b=j; }
            }
        }
        for (int i = 0; i < GAS_PARTICLES; ++i) {
            particles[i].x += particles[i].vx*soonest;
            particles[i].y += particles[i].vy*soonest;
        }
        remaining -= soonest;
        if (hit_a < 0) break;
        Particle *p = &particles[hit_a];
        if (hit_b >= 0) {
            Particle *q = &particles[hit_b];
            float dx=q->x-p->x, dy=q->y-p->y;
            float d2=dx*dx+dy*dy;
            float projection=((q->vx-p->vx)*dx+(q->vy-p->vy)*dy)/d2;
            /* Equal masses exchange normal velocity; tangent stays unchanged. */
            p->vx += projection*dx; p->vy += projection*dy;
            q->vx -= projection*dx; q->vy -= projection*dy;
            ++gas_diagnostics.pair_collisions;
        } else {
            if (wall_axis == 0) {
                p->x = p->vx > 0 ? XMAX : XMIN;
                impulse_x += 2*GAS_MASS*fabsf(p->vx); p->vx = -p->vx;
            } else {
                p->y = p->vy > 0 ? YMAX : YMIN;
                impulse_y += 2*GAS_MASS*fabsf(p->vy); p->vy = -p->vy;
            }
            ++gas_diagnostics.wall_collisions;
        }
    }
    float k = energy();
    gas_diagnostics.kinetic_energy = k;
    gas_diagnostics.temperature = k/(GAS_PARTICLES*GAS_BOLTZMANN);
    pressure_time += GAS_DT; energy_integral += k*GAS_DT;
    if (pressure_time >= GAS_PRESSURE_SECONDS) {
        float px=impulse_x/(2*HEIGHT*pressure_time);
        float py=impulse_y/(2*WIDTH*pressure_time);
        gas_diagnostics.pressure_x=px; gas_diagnostics.pressure_y=py;
        gas_diagnostics.pressure=.5f*(px+py);
        gas_diagnostics.ideal_pressure=energy_integral/(pressure_time*WIDTH*HEIGHT);
        gas_diagnostics.compressibility=gas_diagnostics.pressure/gas_diagnostics.ideal_pressure;
        impulse_x=impulse_y=pressure_time=energy_integral=0;
    }
    ++gas_diagnostics.steps;
}

#ifndef GAS_HOST_TEST
typedef struct { short x, y; } Position;
static Position snapshot[GAS_PARTICLES];
static mutex_t snapshot_mutex;
static void publish_positions(void) {
    Position next[GAS_PARTICLES];
    for (int i=0; i<GAS_PARTICLES; ++i) {
        next[i].x=(short)(particles[i].x+.5f);
        next[i].y=(short)(particles[i].y+.5f);
    }
    mutex_enter_blocking(&snapshot_mutex);
    memcpy(snapshot,next,sizeof snapshot);
    mutex_exit(&snapshot_mutex);
}
static void core1_render(void) {
    Position frame[GAS_PARTICLES];
    while (true) {
        /* Copy before waiting for swap so the rendering window never waits on
         * physics. Only this core touches the VGA drawing API. */
        mutex_enter_blocking(&snapshot_mutex);
        memcpy(frame,snapshot,sizeof frame);
        mutex_exit(&snapshot_mutex);
        while (!draw_start_signal()) tight_loop_contents();
        uint32_t start=time_us_32();
        clearLowFrame(0,BLACK);
        // Redraw text in each back buffer, using Lab 2's 16x32 and 6x7 fonts.
        // Wrap the subtitle to keep it entirely inside the left column.
        drawTextGrotesk32(24,198,"OUTBREAK",WHITE,BLACK);
        drawTextAscii(24,244,"ECE 5730",WHITE,BLACK);
        drawTextAscii(24,260,"Distributed Computing Project",WHITE,BLACK);
        drawRect(GAS_LEFT,GAS_TOP,GAS_RIGHT-GAS_LEFT+1,GAS_BOTTOM-GAS_TOP+1,WHITE);
        for (int i=0; i<GAS_PARTICLES; ++i)
            drawCircle(frame[i].x,frame[i].y,GAS_RADIUS,WHITE);
        uint32_t elapsed=time_us_32()-start;
        if (elapsed>gas_diagnostics.max_render_us) gas_diagnostics.max_render_us=elapsed;
    }
}
int main(void) {
    /* Copied Lab 2 PIO divisors are adapted for this stock clock. */
    set_sys_clock_khz(150000,true);
    initVGA();
    mutex_init(&snapshot_mutex);
    rng_state=get_rand_32();
    if (!rng_state) rng_state=0x57301234u;
    if (!initialize_gas()) panic("Gas configuration cannot fit without overlaps");
    publish_positions();
    multicore_launch_core1(core1_render);
    uint64_t deadline=time_us_64();
    while (true) {
        physics_step();
        publish_positions();
        deadline += GAS_STEP_US;
        if (time_us_64() < deadline) sleep_until(from_us_since_boot(deadline));
        else {
            ++gas_diagnostics.late_steps;
            /* Keep fixed physics dt under overload; slow simulated time instead
             * of unbounded catch-up work. Rendering remains independent. */
            deadline=time_us_64();
        }
    }
}
#endif
