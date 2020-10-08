#include <QCoreApplication>
#include <QDebug>
#include "../../libs/fpc1011f3/fpc1011f3.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    int ret;
    int i, j;

    qDebug("Spi_demo");

    // Buffer geral
    uint8_t buf[16];
    // Buffer de captura da imagem
    uint8_t img_buf[RD_FULL_IMAGE_SIZE_RX];
    // Buffer (x,y) final - imagem
    volatile uint8_t img_array[N_LIN][N_COL];

    // Inicializa buffers de imagem
    memset(img_buf, 0, sizeof(img_buf));
    for (j = 0; j < N_LIN; j++) {
        for (i = 0; i < N_COL; i++) {
            img_array[j][i] = 0;
        }
    }

    // Inicializa fpc1011f3
    fpc1011f3_init();

    // Lê os 16 registros do fpc1011f3
    fpc1011f3_le_regs(buf);

    // Imprime os regs do fpc1011f3
    for (uint8_t i = 0; i < 16; i++) {
        printf("%x ", buf[i]);
    }
    puts("");

    // Lê reg de status
     ret = fpc1011f3_le_status();
     if (ret) {
         printf("Ha dados disponiveis!\r\n");
     }

    // Captura a imagem completa (152 x 200 pixels = 30400 bytes)
    fpc1011f3_le_imagem(img_buf);

    // Transfere imagem para array (x,y)
    for (j = 0; j < N_LIN; j++) {
        for (i = 0; i < N_COL; i++) {
            img_array[j][i] = img_buf[i + (N_COL * j)];
        }
    }

    return a.exec();
}
