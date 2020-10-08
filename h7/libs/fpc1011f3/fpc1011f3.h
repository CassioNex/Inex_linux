/*
* fpc1011f3.h
*/

#ifndef FPC1011F3_H
#define FPC1011F3_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/spi/spidev.h>

#include "../common/common_def.h"

//! [definição de constantes]
//! // Parametros
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
// Matriz de pixels - sizes
#define N_LIN    200
#define N_COL    152
//! [definição de constantes]

//! [definição de estruturas públicas do módulo]
//! [definição de estruturas públicas do módulo]

//! [prototipagem de funções públicas]
PUBLIC int fpc1011f3_init(void);
PUBLIC void fpc1011f3_exit(void);
PUBLIC void fpc1011f3_le_regs(u_int8_t *buf);
PUBLIC int fpc1011f3_le_status(void);
PUBLIC void fpc1011f3_le_imagem(u_int8_t *buf);
//! [prototipagem de funções públicas]

//! [Variáveis globais públicas em memória de dados - usar diretiva extern]
extern PUBLIC union unsigned_char fpc1011_flags;
//! [Variáveis globais públicas em memória de dados - usar diretiva extern]

//! [flags]
/*
* bit0: F_FPC1011_TEMP
* bit1:
* bit2:
* bit3:
* bit4:
* bit5:
* bit6:
* bit7:
*/
#define F_FPC1011_TEMP  	fpc1011_flags.bit0
//! [flags]

//! [Variáveis globais públicas em memória de programa - usar diretiva extern]
//! [Variáveis globais públicas em memória de programa - usar diretiva extern]

//! [macros]
//! [macros]

#ifdef __cplusplus
}
#endif

#endif // FPC1011F3_H
