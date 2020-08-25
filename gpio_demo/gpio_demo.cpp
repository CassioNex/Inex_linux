#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <fcntl.h>
#include <signal.h>
#include <pthread.h>
#include <linux/input.h>

#define LED_B_BRIGHTNESS    "/sys/class/leds/led_b/brightness"
#define LED_Y_BRIGHTNESS    "/sys/class/leds/led_y/brightness"
#define INPUT_EVENT         "/dev/input/event2"
#define LED_MAX_SPEED       10
#define PERIOD_COEFF        16000

unsigned int led_speed;
pthread_mutex_t lock;

// Blink LED
static void *LEDMod(void *dummy)
{
    unsigned int led_period;
    int tmp;
    tmp = open(LED_B_BRIGHTNESS, O_WRONLY);
    if (tmp < 0)
        exit(1);
    while (1) {
        pthread_mutex_lock(&lock);
        led_period = (LED_MAX_SPEED - led_speed) * PERIOD_COEFF;
        pthread_mutex_unlock(&lock);

        write(tmp, "1", 2);
        usleep(led_period);
        write(tmp, "0", 2);
        usleep(led_period);
    }
}

int main()
{
    pthread_t pth;
    struct input_event ev;
    int tmp;
    int key_code;
    int size = sizeof(ev);

    // Configure LED
    led_speed = 5;

    // Create thread
    pthread_mutex_init(&lock, NULL);
    pthread_create(&pth, NULL, LEDMod, (void *)"Blinking LED...");

    // Read event2 (buttons)
    tmp = open(INPUT_EVENT, O_RDONLY);
    if (tmp < 0) {
        printf("\nOpen " INPUT_EVENT " failed!\n");
        return 1;
    }
    /* Read and parse event, update global variable */
    while (1) {
        if (read(tmp, &ev, size) < size) {
            printf("\nReading from " INPUT_EVENT " failed!\n");
            return 1;
        }

        // Checa pressionamento das keys
        if (ev.value == 1 && ev.type == 1) {
            key_code = ev.code;
            // lower speed (bt_1)
            if (key_code == KEY_LINEFEED) {
                // Protect from concurrent read/write
                pthread_mutex_lock(&lock);
                if (led_speed > 0)
                    led_speed -= 1;
                pthread_mutex_unlock(&lock);
            // raise speed
            } else if (key_code == KEY_HOME) {
                pthread_mutex_lock(&lock);
                if (led_speed < 9)
                    led_speed += 1;
                pthread_mutex_unlock(&lock);
            }
            printf("Speed: %i\n", led_speed);
            usleep(1000);
        }
    }
}
