TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    foundfeasible.cpp \
    checkedcolumnstree.cpp

HEADERS += \
    vector.h \
    foundfeasible.h \
    checkedcolumnstree.h

unix:!macx: LIBS += -lblas -llapack
