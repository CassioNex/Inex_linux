#include <QCoreApplication>
#include <QDebug>

#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/types.h>

#include "spidev.h"

// Tamanho da imagem em bytes
#define img_size (30400) /*(4096 - 2)*/

// Nome do dispositivo de driver SPI no diretorio /dev
static const char *device = "/dev/spidev0.0";
// Modo de operacao do controlador SPI do SoC
static uint32_t mode;
// Quantidade de bits da palavra de transferencia SPI
static uint8_t bits = 8;
// Velocidade de comunicacao
static uint32_t speed = (32 * 1000 * 1000);
// Delay entre bytes transferidos, caso necessario
static uint16_t delay;

// Funcao que imprime mensagem de erro e aborta o programa
static void pabort(const char *s)
{
    perror(s);
    abort();
}

/*
 * Funcao que realiza uma transferencia full duplex,
 * ou seja, que escreve e ler ao mesmo tempo
 *
 * Parametros ----
 * fd: inteiro retornado pela funcao open()
 * tx: para de escrita para mensagem SPI
 * rx: buffer de leitura de mensagem SPI
 * len: n. de bytes a serem escritos/lidos
 */
static void transfer(int fd, uint8_t *tx, uint8_t *rx, uint len)
{
    int ret;

    // Estrutura que contem as informacoes para a transmissao da mensagem
    struct spi_ioc_transfer tr = {
        .tx_buf = (unsigned long)tx,
        .rx_buf = (unsigned long)rx,
        .len = len,
        .speed_hz = speed,
        .delay_usecs = delay,
        .bits_per_word = bits,
    };

    // Funcao que faz a transferencia full duplex
    ret = ioctl(fd, SPI_IOC_MESSAGE(1), &tr);
    if (ret < 1) {
        pabort("Nao foi possivel enviar mensagem SPI");
    }
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    qDebug("Spi_demo");

    int ret = 0;
    int fd;
    // Buffers de leitura e escrita básica na spi
    uint8_t tx_buf[20], rx_buf[20];
    // Buffer de leitura de imagem
    uint8_t tx_img_buf[img_size+2], rx_img_buf[img_size+2];

    // Abre o dispositivo de driver SPI e retorna um inteiro
    fd = open(device, O_RDWR);
    if (fd < 0) {
        pabort("Erro ao abrir o dispositivo");
    }

    // Escolha outros modos de operacao que o SPI da sua placa suporta
    mode = /*SPI_CPHA | SPI_CPOL |*/ SPI_MODE_0;

    // Configura o modo de operacao
    ret = ioctl(fd, SPI_IOC_WR_MODE32, &mode);
    if (ret == -1) {
        pabort("Erro ao setar o modo do SPI");
    }

    ret = ioctl(fd, SPI_IOC_RD_MODE32, &mode);
    if (ret == -1) {
        pabort("Erro ao setar o modo do SPI");
    }

    // Configura o tamanho da palavra de transferencia SPI para escrita
    ret = ioctl(fd, SPI_IOC_WR_BITS_PER_WORD, &bits);
    if (ret == -1) {
        pabort("Erro ao setar os bits por palavra");
    }

    // Configura o tamanho da palavra de transferencia SPI para leitura
    ret = ioctl(fd, SPI_IOC_RD_BITS_PER_WORD, &bits);
    if (ret == -1) {
        pabort("Erro ao ler os bits por palavra");
    }

    // Configura a maxima velocidade de transferencia para escrita
    ret = ioctl(fd, SPI_IOC_WR_MAX_SPEED_HZ, &speed);
    if (ret == -1) {
        pabort("Erro ao setar a velocidade maxima em HZ");
    }

    // Configura a maxima velocidade de transferencia para leitura
    ret = ioctl(fd, SPI_IOC_RD_MAX_SPEED_HZ, &speed);
    if (ret == -1) {
        pabort("Erro ao ler a velocidade maxima em HZ");
    }

    // Imprime informacoes de configuracao do SPI
    printf("Modo SPI: 0x%x\n", mode);
    printf("bits por palavra: %d\n", bits);
    printf("Maxima velocidade: %d Hz (%d KHz)\n", speed, speed/1000);

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
    transfer(fd, tx_buf, rx_buf, 20);

    // Imprime os regs do fpc1011f3
    for (uint8_t i = 0; i < 16; i++) {
        printf("%x ", rx_buf[i+4]);
    }
    puts("");

    /****************************************************************/
   // Faz o setup do drivc = ffh
    tx_buf[0] = 0x75;
    tx_buf[1] = 0xff;
    transfer(fd, tx_buf, rx_buf, 2);

    /****************************************************************/
    // Faz a transferencia full duplex - leitura dos regs do fpc1011f3
    tx_buf[0] = 0x50;
    tx_buf[1] = 0x00;
    tx_buf[2] = 0x20;
    tx_buf[3] = 0x00;
    transfer(fd, tx_buf, rx_buf, 20);

    // Verifica novamente os regs do fpc1011f3
    for (uint8_t i = 0; i < 16; i++) {
        printf("%x ", rx_buf[i+4]);
    }
    puts("");

    /****************************************************************/
    // Lê status do fpc1011f3
    tx_buf[0] = 0x21;
    tx_buf[1] = 0x00;
    tx_buf[2] = 0x00;
    transfer(fd, tx_buf, rx_buf, 3);

    // Checa status
    if (rx_buf[2] & 0x01) {
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
    transfer(fd, tx_buf, rx_buf, 2);
    // Fase 2 - Checa status até que dados estejam disponiveis
    tx_buf[0] = 0x21;
    tx_buf[1] = 0x00;
    tx_buf[2] = 0x00;
    transfer(fd, tx_buf, rx_buf, 3);
    while (!(rx_buf[2] & 0x01)) {
        transfer(fd, tx_buf, rx_buf, 3);
    };
    // Fase 3 - Leitura da imagem
    tx_img_buf[0] = 0x20;
    transfer(fd, tx_img_buf, rx_img_buf, (img_size+2));

    // Fecha o dispositivo de driver SPI
    close(fd);

    return a.exec();
}
