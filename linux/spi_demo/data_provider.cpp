#include <QtCore/QFile>
#include <QDebug>
#include "data_provider.h"

// Construtora: faz a conexão signal-slot e dispara o timer
DataProvider::DataProvider(QObject *parent) : QObject(parent)
{
    // Inicializa buffer de imagem
    memset(img_buf, 0, sizeof(img_buf));

    // Inicializa fpc1011f3
    fpc1011f3_init();

    QObject::connect(&timer, &QTimer::timeout, this, &DataProvider::handleTimer);

    timer.setInterval(1000);
    timer.start();
}

// Destrutora
DataProvider::~DataProvider()
{
    fpc1011f3_exit();
}

void DataProvider::handleTimer()
{
    // Debug vars
    // int ret;
    // Buffer geral
    uint8_t buf[16];

    /* Debug
    // Inicializa buffer de imagem
    memset(img_buf, 0, sizeof(img_buf));
    // Inicializa fpc1011f3
    fpc1011f3_init();
    */

    // Debug - Lê os 16 registros do fpc1011f3
    // TODO - Em teoria, não precisaria fazer isso aqui,
    // mas sem esse "gasto" de tempo, o sensor não
    // ativa o status de dados disponíveis, e a
    // fpc1011f3_le_imagem() trava aguardando.
    fpc1011f3_le_regs(buf);

    /* Debug - Imprime os regs do fpc1011f3
    for (uint8_t i = 0; i < 16; i++) {
        printf("%x ", buf[i]);
    }
    puts("");
    */

    /* Debug - Lê reg de status
     ret = fpc1011f3_le_status();
     if (ret) {
         printf("Ha dados disponiveis!\r\n");
     }
     */

    // Captura a imagem completa (152 x 200 pixels = 30400 bytes)
    fpc1011f3_le_imagem(img_buf);

    /* Debug
    qDebug() << "Img captured!";
    */

    emit img_captured(img_buf);
}
