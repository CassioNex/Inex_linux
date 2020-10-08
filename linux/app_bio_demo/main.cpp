#include <QApplication>
#include <QDebug>
#include "window.h"
#include "data_provider.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    DataProvider dp;
    Window window;

    QObject::connect(&dp, &DataProvider::img_captured, &window, &Window::handleImgCaptured);

    window.setFixedSize(480, 800);
    window.setStyleSheet("background-color: white;");
    window.show();

    return app.exec();
}
