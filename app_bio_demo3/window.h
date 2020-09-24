#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QCheckBox>
#include <QPushButton>
#include <QtCore/QElapsedTimer>
#include <unistd.h>
#include "data_provider.h"
#include "finger_img.h"

// Forward declaration
class DataProvider;
class FingerImg;

class Window : public QWidget
{
    Q_OBJECT

public:
    explicit Window(QWidget *parent = nullptr);

signals:

public slots:
    void handleImgCaptured(uint8_t *buf);

private slots:
    void cb_mode_clicked();
    void bt_enroll_ButtonClicked();
    void bt_search_ButtonClicked();

private:
     DataProvider *dp;
     FingerImg  *finger_img;
     QCheckBox  *cb_mode;
     QPushButton *bt_enroll;
     QPushButton *bt_search;
     QElapsedTimer *timer;
     QLabel *vl_enroll_cnt;
};

#endif // WINDOW_H
