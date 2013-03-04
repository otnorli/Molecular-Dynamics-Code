TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    atom.cpp \
    cell.cpp \
    generatequantities.cpp \
    potentials.cpp \
    lib.cpp

HEADERS += \
    atom.h \
    cell.h \
    generatequantities.h \
    potentials.h \
    lib.h

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

