#include <QCoreApplication>
#include <QDebug>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/spi/spidev.h>

// Parametros
// Dispositivo spidev
#define DEV_SPI                 "/dev/spidev0.0"
// Tamanho da imagem em bytes
#define IMG_SIZE                (30400)
// Tamanho da mensagem spi para leitura dos 16 regs
#define RD_ALL_REGS_SIZE_TX     (4)
#define RD_ALL_REGS_SIZE_RX     (16)
// Tamanho da mensagem spi para setup do drivc reg
#define WR_DRIVC_SIZE_TX        (2)
#define WR_DRIVC_SIZE_RX        (0)
// Tamanho da mensagem spi para leitura de status reg
#define RD_STATUS_SIZE_TX       (2)
#define RD_STATUS_SIZE_RX       (1)
// Tamanho da mensagem spi para disparo de read_sensor cmd
#define WR_RDSENSOR_SIZE_TX     (2)
#define WR_RDSENSOR_SIZE_RX     (0)
// Tamanho da mensagem spi para leitura de full image
#define RD_FULL_IMAGE_SIZE_TX   (2)
#define RD_FULL_IMAGE_SIZE_RX   (IMG_SIZE)
// Modo de operação
#define SPI_MODE                SPI_MODE_0
// Velocidade
#define SPI_SPEED               (32*1000*1000)
// N. de bits por palavra
#define SPI_N_BITS              8

// Funcao que realiza uma transferencia spi full duplex
// Parametros:
// fd: file descriptor retornado pela funcao open()
// tx: buffer de escrita
// rx: buffer de leitura
// len_tx: n. de bytes a serem escritos
// len_rx: n. de bytes a serem lidos
static void transfer(int fd, uint8_t *tx, uint8_t *rx, uint len_tx, uint len_rx)
{
    struct spi_ioc_transfer xfer[2] = {0,0};

    xfer[0].tx_buf = (unsigned long)tx;
    xfer[0].len = len_tx;

    xfer[1].tx_buf = 0;
    xfer[1].rx_buf = (unsigned long)rx;
    xfer[1].len = len_rx;

    if (len_tx) {
        if (len_rx) {
            if (ioctl(fd, SPI_IOC_MESSAGE(2), xfer) < 0) {
                perror("SPI_IOC_MESSAGE");
            }
        }
        else {
            if (ioctl(fd, SPI_IOC_MESSAGE(1), xfer) < 0) {
                perror("SPI_IOC_MESSAGE");
            }
        }
    }
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    qDebug("Spi_demo");

    int ret;
    int fd;
    unsigned int mode, speed;
    uint8_t n_bits;
    // Buffers de leitura e escrita básica na spi
    uint8_t tx_buf[RD_ALL_REGS_SIZE_TX], rx_buf[RD_ALL_REGS_SIZE_RX];
    // Buffer de leitura de imagem
    uint8_t tx_img_buf[RD_FULL_IMAGE_SIZE_TX], rx_img_buf[RD_FULL_IMAGE_SIZE_RX];

    // open device node
    fd = open(DEV_SPI, O_RDWR);
    if (fd < 0) {
        printf("ERROR open %s ret=%d\n", DEV_SPI, fd);
        return -1;
    }

    // set spi mode
    mode = SPI_MODE;
    if (ioctl(fd, SPI_IOC_WR_MODE32, &mode) < 0) {
        printf("ERROR ioctl() set mode\n");
        return -1;
    }
    if (ioctl(fd, SPI_IOC_RD_MODE32, &ret) < 0) {
        printf("ERROR ioctl() get mode\n");
        return -1;
    }
    else {
        printf("mode set to %d\n", (unsigned int)ret);
    }

    // Configura o tamanho da palavra de transferencia SPI
    n_bits = SPI_N_BITS;
    if (ioctl(fd, SPI_IOC_WR_BITS_PER_WORD, &n_bits) < 0) {
        printf("ERROR ioctl() set n_bits\n");
        return -1;
    }
    if (ioctl(fd, SPI_IOC_RD_BITS_PER_WORD, &ret) < 0) {
        printf("ERROR ioctl() get n_bits\n");
        return -1;
    }
    else {
        printf("n_bits set to %d\n", (unsigned int)ret);
    }

    // set spi speed
    speed = SPI_SPEED;
    if (ioctl(fd, SPI_IOC_WR_MAX_SPEED_HZ, &speed) < 0) {
        printf("ERROR ioctl() set speed\n");
        return -1;
    }
    if (ioctl(fd, SPI_IOC_RD_MAX_SPEED_HZ, &ret) < 0) {
        printf("ERROR ioctl() get speed\n");
        return -1;
    }
    else {
        printf("speed set to %d\n", ret);
    }

    // Inicializa buffers de tx e rx
    memset(tx_buf, 0, sizeof(tx_buf));
    memset(rx_buf, 0, sizeof(rx_buf));
    memset(tx_img_buf, 0, sizeof(tx_img_buf));
    memset(rx_img_buf, 0, sizeof(rx_img_buf));

    /****************************************************************/
    // Faz a transferencia full duplex - leitura dos regs do fpc1011f3
    tx_buf[0] = 0x50;
    tx_buf[1] = 0x00;
    tx_buf[2] = 0x20;
    tx_buf[3] = 0x00;
    transfer(fd, tx_buf, rx_buf, RD_ALL_REGS_SIZE_TX, RD_ALL_REGS_SIZE_RX);

    // Imprime os regs do fpc1011f3
    for (uint8_t i = 0; i < 16; i++) {
        printf("%x ", rx_buf[i]);
    }
    puts("");

    /****************************************************************/
   // Faz o setup do drivc = ffh
    tx_buf[0] = 0x75;
    tx_buf[1] = 0xff;
    transfer(fd, tx_buf, rx_buf, WR_DRIVC_SIZE_TX, WR_DRIVC_SIZE_RX);

    /****************************************************************/
    // Faz a transferencia full duplex - leitura dos regs do fpc1011f3
    tx_buf[0] = 0x50;
    tx_buf[1] = 0x00;
    tx_buf[2] = 0x20;
    tx_buf[3] = 0x00;
    transfer(fd, tx_buf, rx_buf, RD_ALL_REGS_SIZE_TX, RD_ALL_REGS_SIZE_RX);

    // Verifica novamente os regs do fpc1011f3
    for (uint8_t i = 0; i < 16; i++) {
        printf("%x ", rx_buf[i]);
    }
    puts("");

    /****************************************************************/
    // Lê status do fpc1011f3
    tx_buf[0] = 0x21;
    tx_buf[1] = 0x00;
    transfer(fd, tx_buf, rx_buf, RD_STATUS_SIZE_TX, RD_STATUS_SIZE_RX);

    // Checa status
    if (rx_buf[0] & 0x01) {
        printf("Ha dados disponiveis!\r\n");
    }
    else {
        printf("Nao Ha dados disponiveis!\r\n");
    }

    /****************************************************************/
    // Leitura da imagem
    // Fase 1 - envia comando rd_sensor
    tx_buf[0] = 0x11;
    tx_buf[1] = 0x00;
    transfer(fd, tx_buf, rx_buf, WR_RDSENSOR_SIZE_TX, WR_RDSENSOR_SIZE_RX);
    // Fase 2 - Checa status até que dados estejam disponiveis
    tx_buf[0] = 0x21;
    tx_buf[1] = 0x00;
    transfer(fd, tx_buf, rx_buf, RD_STATUS_SIZE_TX, RD_STATUS_SIZE_RX);
    while (!(rx_buf[0] & 0x01)) {
        transfer(fd, tx_buf, rx_buf, RD_STATUS_SIZE_TX, RD_STATUS_SIZE_RX);
    };
    // Fase 3 - Leitura da imagem
    tx_img_buf[0] = 0x20;
    tx_img_buf[1] = 0x00;
    transfer(fd, tx_img_buf, rx_img_buf, RD_FULL_IMAGE_SIZE_TX, RD_FULL_IMAGE_SIZE_RX);

    // Fecha o dispositivo de driver SPI
    close(fd);

    return a.exec();
}
