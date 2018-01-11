TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    smartCalculation.cpp \
    multizoneFit.cpp \
    mobilityspectrum.cpp \
    MagneticFieldDependence.cpp \
    FilterParams.cpp \
    FilteringUnit.cpp \
    ExtrapolateUnit.cpp \
    DataSaver.cpp \
    commonFunctions.cpp \
    ErrorEst/NoiseUnit.cpp \
    ErrorEst/film.cpp \
    ErrorEst/clMagneticFieldDependences.cpp

DISTFILES += \
    MCTConsole.pro.user \
    MCTParamsExplorer.sublime-workspace \
    MCTParamsExplorer.sublime-project

HEADERS += \
    smartCalculation.h \
    multizoneFit.h \
    mobilityspectrum.h \
    main.h \
    MagneticFieldDependence.h \
    FilterParams.h \
    FilteringUnit.h \
    ExtrapolateUnit.h \
    DataSaver.h \
    commonFunctions.h \
    ErrorEst/NoiseUnit.h \
    ErrorEst/film.h \
    ErrorEst/clMagneticFieldDependences.h
