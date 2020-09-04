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

    // Define os layouts
    QHBoxLayout *hlay   = new QHBoxLayout;
    QVBoxLayout *layout = new QVBoxLayout;

    // Define as fontes de texto
    QFont ft;
    ft.setPointSize(10);
    ft.setBold(false);

    QFont fv;
    fv.setPointSize(10);
    fv.setBold(true);

    QLabel *lb_score = new QLabel(tr("Score: "));
    lb_score->setAlignment(Qt::AlignHCenter);
    lb_score->setFont(ft);

    vl_score = new QLabel;
    vl_score->setAlignment(Qt::AlignHCenter);
    vl_score->setFont(fv);

    QLabel *lb_index = new QLabel(tr("Pos: "));
    lb_index->setAlignment(Qt::AlignHCenter);
    lb_index->setFont(ft);

    vl_index = new QLabel;
    vl_index->setAlignment(Qt::AlignHCenter);
    vl_index->setFont(fv);

    // Adiciona os widgets aos layouts
    hlay->addWidget(lb_score);
    hlay->addWidget(vl_score);
    hlay->addWidget(lb_index);
    hlay->addWidget(vl_index);

    // Objeto para exibir imagem do fingerprint
    finger = new QLabel();
    finger->setAlignment(Qt::AlignHCenter);
    finger->setPixmap(QPixmap(N_COL, N_LIN));

    layout->addLayout(hlay);
    layout->addWidget(finger);

    setLayout(layout);

    // Inicializa flag de enroll x search
    m_mode = ENROLL_MODE;

    // Abre descriptors para os leds de enroll e match
    m_fd_led_ok = open(LED_B_BRIGHTNESS, O_WRONLY);
    if (m_fd_led_ok < 0) {
        exit(1);
    }
    m_fd_led_error = open(LED_Y_BRIGHTNESS, O_WRONLY);
    if (m_fd_led_error < 0) {
        exit(1);
    }
    // Default - leds apagados
    write(m_fd_led_ok, "0", 2);
    write(m_fd_led_error, "0", 2);

    // Inicializa contador de templates cadastrados
    enroll_cnt = 0;
}

void FingerImg::handleImgCaptured(uint8_t *buf)
{
    /* Debug
    uint8_t img_buf[IMG_SIZE];
    */
    int search_res;
    int match_res;
    static int match_res_n1;


    /* Não aloca no stack, usa heap (ver construtora)
    static uint8_t feature_enroll[MAX_FEATUREVECT_LEN];
    static uint8_t feature_search[MAX_FEATUREVECT_LEN];
    static uint8_t feature_db[(N_TEMPLATES + 1) * MAX_FEATUREVECT_LEN];
    static uint8_t *p_db = feature_db;
    */
    // Alocação de memória dinâmica para os templates
    static uint8_t *pfeature_enroll;
    static uint8_t *pfeature_search;
    static uint8_t *pfeature_db;
    static uint8_t *pdb_enroll;
    static uint8_t *pdb_search;
    static bool flag = false;

    if (!flag) {
        pfeature_enroll = (uint8_t *)malloc(sizeof(uint8_t)*MAX_FEATUREVECT_LEN);
        pfeature_search = (uint8_t *)malloc(sizeof(uint8_t)*MAX_FEATUREVECT_LEN);
        pfeature_db     = (uint8_t *)malloc(sizeof(uint8_t)*MAX_FEATUREVECT_LEN*N_TEMPLATES);
        pdb_enroll      = pfeature_db;
        pdb_search      = pfeature_db;
        memset(pfeature_db, 0, sizeof(uint8_t)*MAX_FEATUREVECT_LEN*N_TEMPLATES);
        flag = true;
    }

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

    // Enroll mode - salva template no database buffer
    if (m_mode == ENROLL_MODE) {
        // Apaga score e index
        vl_score->setText(NULL);
        vl_index->setText(NULL);
        // Extrai o template (feature extraction)
        if (create_template((BYTE*)/*img_buf*/buf, N_COL, N_LIN, (BYTE*)pfeature_enroll) == 1) {
           qDebug() << "Template extracted!";
           // Salva template no database buffer
           memcpy(pdb_enroll, pfeature_enroll, MAX_FEATUREVECT_LEN);
           pdb_enroll += MAX_FEATUREVECT_LEN;
           // Incrementa contador de templates aprendidos
           enroll_cnt++;
        }
        else {
            qDebug() << "Error on template extraction!";
        }
    }
    // Search mode
    else if (m_mode == SEARCH_MODE) {
        // Extrai o template (feature extraction)
        if (create_template((BYTE*)/*img_buf*/buf, N_COL, N_LIN, (BYTE*)pfeature_search) == 1) {
           qDebug() << "Template extracted!";
           // Faz a busca no database
#ifdef USA_SEARCH
    #if (VERSAO_ALG != ALG_V2_0)
           search_res = finger_search((BYTE *)pfeature_search, (BYTE *)pfeature_db, N_TEMPLATES, MEDIUM_LEVEL);
    #else
           search_res = finger_search((BYTE *)pfeature_search, (BYTE *)pfeature_db, N_TEMPLATES, MEDIUM_LEVEL, (BYTE *)/*img_buf*/buf);
    #endif
#else
           // Usando apenas finger_match
           // Reposiciona pointer no início do db
           pdb_search = pfeature_db;
           // Inicializa search_res
           search_res = -1;
           // Inicializa match_res anterior
           match_res_n1 = 0;
           for (int i = 0; i < N_TEMPLATES; i++) {
    #if (VERSAO_ALG != ALG_V2_0)
                match_res = finger_match((BYTE *)pfeature_search, (BYTE *)pdb_search, MEDIUM_LEVEL);
    #else
                match_res = finger_match((BYTE *)pfeature_search, (BYTE *)pdb_search, MEDIUM_LEVEL, (BYTE *)/*img_buf*/buf);
    #endif

    #if (VERSAO_ALG == ALG_V_PC)
                if (match_res > match_res_n1) {
                    match_res_n1 = match_res;
                    if (match_res >= TH_MATCH) {
                        search_res = i;
                    }
                }
                // Próxima busca
                pdb_search += MAX_FEATUREVECT_LEN;
           }
    #else
                if (match_res == 1) {
                    search_res = i;
                    break;
                }
                // Próxima busca
                pdb_search += MAX_FEATUREVECT_LEN;
           }
    #endif
#endif
           if (search_res < 0) {
               qDebug() << "Fingerprint not in database!";
               write(m_fd_led_ok, "0", 2);
               write(m_fd_led_error, "1", 2);
           }
           else {
               qDebug() << "Fingerprint in database position: " << search_res << "!";
               write(m_fd_led_ok, "1", 2);
               write(m_fd_led_error, "0", 2);
           }
           // Exibe o score e o index
           vl_score->setText(QString::number(match_res_n1));
           vl_index->setText(QString::number(search_res));
        }
        else {
            qDebug() << "Error on template extraction!";
            vl_score->setText(NULL);
            vl_index->setText(NULL);
        }
    }
}
