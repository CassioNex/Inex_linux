#include <QtWidgets>
#include "finger_img.h"

FingerImg::FingerImg(QWidget *parent) : QWidget(parent)
{
    /* Debug
    int i, j;

    // Inicializa buffer de imagem
    for (j = 0; j < N_LIN; j++) {
        for (i = 0; i < N_COL; i++) {
            img_array[j][i] = 0;
        }
    }
    */

    QVBoxLayout *layout = new QVBoxLayout;

    QLabel *title = new QLabel(tr("Finger Image"));
    title->setAlignment(Qt::AlignHCenter);

    finger = new QLabel();
    finger->setAlignment(Qt::AlignHCenter);
    finger->setPixmap(QPixmap(N_COL, N_LIN));

    layout->addWidget(title);
    layout->addWidget(finger);

    setLayout(layout);

    // Inicializa flag de enrol x match
    m_enroll_match = ENROL_MODE;
}

void FingerImg::handleImgCaptured(uint8_t *buf)
{
    /* Debug
    uint8_t img_buf[IMG_SIZE];
    */
    static uint8_t feature_enrol[MAX_FEATUREVECT_LEN];
    static uint8_t feature_match[MAX_FEATUREVECT_LEN];
    int match_res;

    /* Debug
    volatile int i, j, k;
    // Transfere imagem para array (x,y)
    for (j = 0; j < N_LIN; j++) {
        for (i = 0; i < N_COL; i++) {
            img_array[j][i] = *buf++;
        }
    }
    */

    /* Debug - Transfere imagem para img_buf
    for (int i = 0; i < IMG_SIZE; i++) {
        img_buf[i] = *buf++;
    }
    */

    // Carrega buffer de imagem em objeto QImage,
    // e plota via QPixmap no label finger
    QImage img(/*img_buf*/ buf, N_COL, N_LIN, QImage::Format_Grayscale8);
    finger->setPixmap(QPixmap::fromImage(img));

    // Enrol mode
    if (m_enroll_match == ENROL_MODE) {
        // Extrai o template (feature extraction)
        if (create_template((BYTE*)/*img_buf*/buf, N_COL, N_LIN, (BYTE*)feature_enrol) == 1) {
           qDebug() << "Template extracted!";
        }
        else {
            qDebug() << "Error on template extraction!";
        }
    }
    // Match mode
    else {
        // Extrai o template (feature extraction)
        if (create_template((BYTE*)/*img_buf*/buf, N_COL, N_LIN, (BYTE*)feature_match) == 1) {
           qDebug() << "Template extracted!";
        }
        else {
            qDebug() << "Error on template extraction!";
        }

        // Checa se o template extraído confere com o enrolled
#if (VERSAO_ALG != ALG_V2_0)
        match_res = finger_match((BYTE *)feature_match, (BYTE *)feature_enrol, LOW_LEVEL);
#else
        match_res = finger_match((BYTE *)feature_match, (BYTE *)feature_enrol, LOW_LEVEL, (BYTE *)/*img_buf*/buf);
#endif
        if (match_res == 1) {
           qDebug() << "Matched!";
        }
        else {
            qDebug() << "Error on match!";
        }
    }
}
