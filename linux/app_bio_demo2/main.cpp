#include <QApplication>
#include <QDebug>
#include "window.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    Window window;

    window.setFixedSize(480, 800);
    window.setStyleSheet("background-color: white;");
    window.show();

    return app.exec();
}
