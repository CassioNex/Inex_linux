#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/spi/spidev.h>

#define DEV_SPI "/dev/spidev0.0"

int main(int argc, char *argv[])
{
    int fd;
    int ret;
    unsigned int mode, speed;
    char tx_buf[1];
    char rx_buf[1];
    struct spi_ioc_transfer xfer[2] = {0};

    // open device node
    fd = open(DEV_SPI, O_RDWR);
    if (fd < 0) {
        printf("ERROR open %s ret=%d\n", DEV_SPI, fd);
        return -1;
    }

    // set spi mode
    mode = SPI_MODE_0;
    if (ioctl(fd, SPI_IOC_WR_MODE32, &mode) < 0) {
        printf("ERROR ioctl() set mode\n");
        return -1;
    }
    if (ioctl(fd, SPI_IOC_RD_MODE32, &ret) < 0) {
        printf("ERROR ioctl() get mode\n");
        return -1;
    } else
        printf("mode set to %d\n", (unsigned int)ret);

    // set spi speed
    speed = 25*1000*1000;
    if (ioctl(fd, SPI_IOC_WR_MAX_SPEED_HZ, &speed) < 0) {
        printf("ERROR ioctl() set speed\n");
        return -1;
    }
    if (ioctl(fd, SPI_IOC_RD_MAX_SPEED_HZ, &ret) < 0) {
        printf("ERROR ioctl() get speed\n");
        return -1;
    } else
        printf("speed set to %d\n", ret);

    // transfer data
    tx_buf[0] = 0xa5;
    xfer[0].tx_buf = (unsigned long)tx_buf;
    xfer[0].len = 1;
    xfer[1].rx_buf = (unsigned long)rx_buf;
    xfer[1].len = 1;

    do {
        if (ioctl(fd, SPI_IOC_MESSAGE(2), xfer) < 0)
            perror("SPI_IOC_MESSAGE");
        usleep(100*1000);
    } while (1);

    // close device node
    close(fd);

    return 0;
}