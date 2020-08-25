#ifndef FINGER_IMG_H
#define FINGER_IMG_H

#include <QWidget>
#include "data_provider.h"
#include "../libs/common/common_def.h"

#if (VERSAO_ALG == ALG_V_PC)
    #include "../libs/finger_alg/v_pc/fingerapi.h"
#elif (VERSAO_ALG == ALG_V1_0)
    #include "../libs/finger_alg/v1_0/fingerapi.h"
#elif (VERSAO_ALG == ALG_V2_0)
    #include "../libs/finger_alg/v2_0/fingerapi.h"
#endif

// Constantes
// Modo
#define ENROLL_MODE 0
#define SEARCH_MODE 1

// N. de templates (database em ram)
#define N_TEMPLATES 9000

// Leds
#define LED_B_BRIGHTNESS    "/sys/class/leds/led_b/brightness"
#define LED_Y_BRIGHTNESS    "/sys/class/leds/led_y/brightness"

// Forward declaration
class QLabel;

class FingerImg : public QWidget
{
    Q_OBJECT

public:
    explicit FingerImg(QWidget *parent = nullptr);
    bool m_mode;
    /* Debug - Buffer (x,y) final - imagem
    uint8_t img_array[N_LIN][N_COL];
    */
    int enroll_cnt;

signals:

public slots:
    void handleImgCaptured(uint8_t *buf);

private:
    QLabel *finger;
    int m_fd_led_ok;
    int m_fd_led_error;
};

#endif // FINGER_IMG_H
