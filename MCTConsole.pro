TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

LIBS += \
       -lboost_system\

SOURCES += main.cpp \
    smartCalculation.cpp \
    multizoneFit.cpp \
    mobilityspectrum.cpp \
  #  mobilityspectrumBoost.cpp \
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
  #  mobilityspectrumBoost.h \
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
