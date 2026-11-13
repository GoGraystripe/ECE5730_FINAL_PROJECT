// Two-Pico UART link (same code on both boards)
//   Serial monitor: USB or uart0 (GP0/GP1)
//   Link: uart1  GP4 (TX) -> other GP5 (RX), GP5 (RX) <- other GP4 (TX), common GND
//   Commands: toggle, led on, led off

#include <stdio.h>
#include <string.h>
#include "pico/stdlib.h"
#include "hardware/uart.h"

#define LINK_UART uart1
#define LED_PIN   PICO_DEFAULT_LED_PIN
#define BUF_SIZE  64

static void handle(const char *cmd)
{
    static bool led = false;

    if      (!strcmp(cmd, "toggle"))  led = !led;
    else if (!strcmp(cmd, "led on"))  led = true;
    else if (!strcmp(cmd, "led off")) led = false;
    else { printf("[remote] unknown: %s\n", cmd); return; }

    gpio_put(LED_PIN, led);
    printf("[remote] %s -> LED %s\n", cmd, led ? "on" : "off");
}

int main()
{
    stdio_init_all();

    uart_init(LINK_UART, 115200);
    gpio_set_function(4, GPIO_FUNC_UART);
    gpio_set_function(5, GPIO_FUNC_UART);

    gpio_init(LED_PIN);
    gpio_set_dir(LED_PIN, GPIO_OUT);

    char in[BUF_SIZE],  rx[BUF_SIZE];
    int  n_in = 0,      n_rx = 0;

    while (true) {
        // typed on serial monitor -> send over link
        int c = getchar_timeout_us(0);
        if (c != PICO_ERROR_TIMEOUT) {
            if (c == '\r' || c == '\n') {
                if (n_in) {
                    in[n_in] = 0;
                    uart_puts(LINK_UART, in);
                    uart_putc_raw(LINK_UART, '\n');
                    printf("\n[sent] %s\n", in);
                    n_in = 0;
                }
            } else if (n_in < BUF_SIZE - 1) {
                in[n_in++] = c;
                putchar(c);  // echo
            }
        }

        // received over link -> run command
        while (uart_is_readable(LINK_UART)) {
            char r = uart_getc(LINK_UART);
            if (r == '\n') {
                rx[n_rx] = 0;
                if (n_rx) handle(rx);
                n_rx = 0;
            } else if (n_rx < BUF_SIZE - 1) {
                rx[n_rx++] = r;
            }
        }
    }
}
