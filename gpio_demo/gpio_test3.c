#include <stdio.h>
#include <unistd.h>
#include <gpiod.h>

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <linux/input.h>
 
int main(int argc, char *argv[])
{
    int fd ;
    int bytes_lidos;
    struct input_event ev;
    size_t size_ie = sizeof(struct input_event);
    char buf[size_ie];

    fd = open("/dev/input/event2", O_RDONLY);
    bytes_lidos = read(fd, buf, size_ie * 1);

    while (1) {
        sleep(1);
        bytes_lidos = read(fd, buf, size_ie * 1);
        if (ev.code == 101) {
            printf("Tecla bt_1: %d \r\n", ev.value);
        }
        if (ev.code == 102) {
            printf("Tecla bt_2: %d \r\n", ev.value);
        }
    }
    
    return 0;
}
 
