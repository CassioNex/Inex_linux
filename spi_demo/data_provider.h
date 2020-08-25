#ifndef DATAPROVIDER_H
#define DATAPROVIDER_H

#include <QObject>
#include <QtCore/QTimer>
#include "../libs/fpc1011f3/fpc1011f3.h"

class DataProvider : public QObject
{
    Q_OBJECT

public:
    explicit DataProvider(QObject *parent = nullptr);
    ~DataProvider();
    // Buffer de captura da imagem
    uint8_t img_buf[RD_FULL_IMAGE_SIZE_RX];

private slots:
    void handleTimer();

signals:
    void img_captured(uint8_t *buf);

private:
    QTimer timer;
};

#endif // DATAPROVIDER_H
