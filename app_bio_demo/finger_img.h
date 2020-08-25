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
#define ENROL_MODE 0
#define MATCH_MODE 1

// Forward declaration
class QLabel;

class FingerImg : public QWidget
{
    Q_OBJECT

public:
    explicit FingerImg(QWidget *parent = nullptr);
    bool m_enroll_match;
    /* Debug - Buffer (x,y) final - imagem
    uint8_t img_array[N_LIN][N_COL];
    */

signals:

public slots:
    void handleImgCaptured(uint8_t *buf);

private:
    QLabel *finger;
};

#endif // FINGER_IMG_H
