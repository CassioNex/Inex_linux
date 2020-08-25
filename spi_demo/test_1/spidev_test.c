/*
 * Exemplo de software para comunicacao SPI
 * 
 * Autor: Vinicius Maciel
 *
 * comando para compilacao: gcc spi_teste.c -o spi_teste
 * 
 * comando para cross compilacao: arm-linux-gnueabihf-gcc spi_teste.c -o spi_teste
 */
 
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/types.h>
 
#include "spidev.h"
 
#define ARRAY_SIZE(a) (sizeof(a) / sizeof((a)[0]))
 
/* Funcao que imprime mensagem de erro e aborta o programa
 * 
 */
static void pabort(const char *s)
{
    perror(s);
    abort();
}
 
// Nome do dispositivo de driver SPI no diretorio /dev 
static const char *device = "/dev/spidev0.0"; // (Mude de acordo com sua placa)
// Modo de operacao do controlador SPI do SoC
static uint32_t mode;
// Quantidade de bits da palavra de transferencia SPI
static uint8_t bits = 8;
// Velocidade de comunicao
static uint32_t speed = 1000000;
// Delay entre bytes transferidos, caso necessario
static uint16_t delay;
 
/*
 * Funcao que realiza uma transferencia full duplex,
 * ou seja, que escreve e ler ao mesmo tempo
 * 
 * Parametros ----
 * fd: inteiro retornado pela funcao open()
 * tx: para de escrita para mensagem SPI
 * rx: buffer de leitura de mensagem SPI
 */
static void transfer(int fd, uint8_t *tx, uint8_t *rx)
{
    int ret;
    
    // Estrutura que contem as informacoes para a transmissao da mensagem
    struct spi_ioc_transfer tr = {
	    .tx_buf = (unsigned long)tx,
	    .rx_buf = (unsigned long)rx,
	    .len = ARRAY_SIZE(tx),
	    .delay_usecs = delay,
	    .speed_hz = speed,
	    .bits_per_word = bits,
    };	
 
    // Funcao que faz a transferencia full duplex
    ret = ioctl(fd, SPI_IOC_MESSAGE(1), &tr);
    if (ret < 1)
	    pabort("Nao foi possivel enviar mensagem SPI");
 
    // Imprimi a mensagem recebida, caso haja alguma
    for (ret = 0; ret < ARRAY_SIZE(rx); ret++) {	    
	    printf("%x ", rx[ret]);
    }
    puts("");
}
 
 
int main(int argc, char *argv[])
{
    int ret = 0, i;
    int fd;
    // Buffers para leitura e escrita da mensagem SPI
    uint8_t tx_buffer[32], rx_buffer[32];
    
 
    // Abri o dispositivo de driver SPI e retorna um inteiro
    fd = open(device, O_RDWR);
    if (fd < 0)
	pabort("Erro ao abrir o dispositivo");
 
    // Escolha outros modos de operacao que o SPI da sua placa suporta
    mode = SPI_CPHA | SPI_CPOL | SPI_MODE_0;
    
    // Configura o modo de operacao
    ret = ioctl(fd, SPI_IOC_WR_MODE32, &mode);
    if (ret == -1)
	    pabort("Erro ao setar o modo do SPI");
 
    ret = ioctl(fd, SPI_IOC_RD_MODE32, &mode);
    if (ret == -1)
	    pabort("Erro ao setar o modo do SPI");
 
    // Configura o tamanho da palavra de transferencia SPI para escrita
    ret = ioctl(fd, SPI_IOC_WR_BITS_PER_WORD, &bits);
    if (ret == -1)
	    pabort("Erro ao setar os bits por palavra");
 
    // Configura o tamanho da palavra de transferencia SPI para leitura
    ret = ioctl(fd, SPI_IOC_RD_BITS_PER_WORD, &bits);
    if (ret == -1)
	    pabort("Erro ao ler os bits por palavra");
 
    // Configura a maxima velocidade de transferencia para escrita
    ret = ioctl(fd, SPI_IOC_WR_MAX_SPEED_HZ, &speed);
    if (ret == -1)
	    pabort("Erro ao setar a velocidade maxima em HZ");
 
    // Configura a maxima velocidade de transferencia para leitura
    ret = ioctl(fd, SPI_IOC_RD_MAX_SPEED_HZ, &speed);
    if (ret == -1)
	    pabort("Erro ao ler a velocidade maxima em HZ");
 
    // Imprimi informacoes de configuracao do SPI
    printf("Modo SPI: 0x%x\n", mode);
    printf("bits por palavra: %d\n", bits);
    printf("Maxima velocidade: %d Hz (%d KHz)\n", speed, speed/1000);
    
    // Preenche o buffer de transferencia com numeros de 0 a 9
    for (i = 0; i < 9; i++)
	tx_buffer[i] = i;
 
    // Faz a transferencia full duplex
    transfer(fd, tx_buffer, rx_buffer);
 
    // Fecha o dispositivo de driver SPI
    close(fd);
 
    return ret;
}