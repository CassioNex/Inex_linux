/*
* fpc1011f3.c
*/

#include "fpc1011f3.h"

//! [prototipagem de funções privadas ao módulo]
PRIVATE void transfer(int fd, u_int8_t *tx, u_int8_t *rx, uint len_tx, uint len_rx);
//! [prototipagem de funções privadas ao módulo]

//! [Variáveis globais públicas em memória de dados]
PUBLIC union unsigned_char fpc1011_flags;
//! [Variáveis globais públicas em memória de dados]

//! [Variáveis globais públicas em memória de programa]
//! [Variáveis globais públicas em memória de programa]

//! [Variáveis globais privadas em memória de dados]
// Spi file descriptor
PRIVATE int fd;
// Spi mode
PRIVATE unsigned int mode;
// SPi speed
PRIVATE unsigned int speed;
// Spi n_bits per word
PRIVATE u_int8_t n_bits;
// Buffers de leitura e escrita básica na spi
PRIVATE u_int8_t tx_buf[RD_ALL_REGS_SIZE_TX];
PRIVATE u_int8_t rx_buf[RD_ALL_REGS_SIZE_RX];
//! [Variáveis globais privadas em memória de dados]

//! [Variáveis globais privadas em memória de programa]
//! [Variáveis globais privadas em memória de programa]

//! [Funções públicas]
/**
 * \external
 * fpc1011f3_init
 *
 * Inicializa spi master e o sensor fpc1011f3
 *
 * \param[in]	void
 * \param[out]	void
 *
 * \return		int (-1: error / 0: ok)
 */
PUBLIC int fpc1011f3_init(void)
{
    int ret;

    // Zera flags do módulo
    fpc1011_flags.value = 0;

    // Inicializa buffers de tx e rx
    memset(tx_buf, 0, sizeof(tx_buf));
    memset(rx_buf, 0, sizeof(rx_buf));

    // Inicializa spi
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

   // Faz o setup do reg drivc = ffh
    tx_buf[0] = 0x75;
    tx_buf[1] = 0xff;
    transfer(fd, tx_buf, rx_buf, WR_DRIVC_SIZE_TX, WR_DRIVC_SIZE_RX);

    return 0;
}

/**
 * \external
 * fpc1011f3_exit
 *
 * Finaliza o acesso à spi e ao sensor
 *
 * \param[in]	void
 * \param[out]	void
 *
 * \return		void
 */
PUBLIC void fpc1011f3_exit(void)
{
    // Fecha o dispositivo de driver SPI
    close(fd);
}

/**
 * \external
 * fpc1011f3_le_regs
 *
 * Lê os registros do fpc1011f3
 *
 * \param[in]	void
 * \param[out]	uint8_t buf
 *
 * \return		void
 */
PUBLIC void fpc1011f3_le_regs(u_int8_t *buf)
{
    tx_buf[0] = 0x50;
    tx_buf[1] = 0x00;
    tx_buf[2] = 0x20;
    tx_buf[3] = 0x00;
    transfer(fd, tx_buf, buf, RD_ALL_REGS_SIZE_TX, RD_ALL_REGS_SIZE_RX);
}

/**
 * \external
 * fpc1011f3_le_status
 *
 * Lê o reg de status do fpc1011f3
 *
 * \param[in]	void
 * \param[out]	void
 *
 * \return		int ret (0: sem dados / 1: com dados disponíveis)
 */
PUBLIC int fpc1011f3_le_status(void)
{
    volatile int ret = 0;
    tx_buf[0] = 0x21;
    tx_buf[1] = 0x00;
    transfer(fd, tx_buf, rx_buf, RD_STATUS_SIZE_TX, RD_STATUS_SIZE_RX);

    // Checa status
    if (rx_buf[0] & 0x01) {
        ret = 1;
    }
    else {
        ret = 0;
    }
    return ret;
}

/**
 * \external
 * fpc1011f3_le_imagem
 *
 * Lê imagem completa do fpc1011f3
 *
 * \param[in]	void
 * \param[out]	uint8_t *buf
 *
 * \return		void
 */
PUBLIC void fpc1011f3_le_imagem(u_int8_t *buf)
{
    // Fase 1 - envia comando rd_sensor
    tx_buf[0] = 0x11;
    tx_buf[1] = 0x00;
    transfer(fd, tx_buf, rx_buf, WR_RDSENSOR_SIZE_TX, WR_RDSENSOR_SIZE_RX);
    // Fase 2 - Checa status até que dados estejam disponiveis
    tx_buf[0] = 0x21;
    tx_buf[1] = 0x00;
    rx_buf[0] = 0x00;
    while (!(rx_buf[0] & 0x01)) {
        transfer(fd, tx_buf, rx_buf, RD_STATUS_SIZE_TX, RD_STATUS_SIZE_RX);
    };
    // Fase 3 - Leitura da imagem
    tx_buf[0] = 0x20;
    tx_buf[1] = 0x00;
    transfer(fd, tx_buf, buf, RD_FULL_IMAGE_SIZE_TX, RD_FULL_IMAGE_SIZE_RX);
}
//! [Funções públicas]

//! [Funções privadas]
/**
 * \internal
 * transfer
 *
 * Funcao que realiza uma transferencia spi full duplex
 * fd: file descriptor retornado pela funcao open()
 * tx: buffer de escrita
 * rx: buffer de leitura
 * len_tx: n. de bytes a serem escritos
 * len_rx: n. de bytes a serem lidos
 *
 * \param[in]	int fd, uint8_t *tx, uint len_tx, uint len_rx
 * \param[out]	uint8_t *rx
 *
 * \return		void
 */
PRIVATE void transfer(int fd, u_int8_t *tx, u_int8_t *rx, uint len_tx, uint len_rx)
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
//! [Funções privadas]
