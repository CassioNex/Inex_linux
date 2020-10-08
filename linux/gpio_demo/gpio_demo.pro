QT += widgets
SOURCES = \
    gpio_demo.cpp
INSTALLS += target
target.path = /usr/bin

unix:!macx: LIBS += -L$$PWD/../buildroot/output/target/usr/lib/ -lgpiod

INCLUDEPATH += $$PWD/../buildroot/output/target/usr
DEPENDPATH += $$PWD/../buildroot/output/target/usr
