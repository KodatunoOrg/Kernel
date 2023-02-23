# Qt 設定 Preferences -> ビルドと実行 -> Default Build Properties -> Default build directory を
# ./build
# に変えておく
CONFIG -= qt
TEMPLATE = lib
CONFIG += staticlib
CONFIG += c++17

# Boost C++ のインクルードパスを各自の環境に合わせて追加
# BoostライブラリはヘッダーのみでOK（ビルドする必要なし）
INCLUDEPATH += \
    ./include/ \
    ../boost

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    src/NURBS.cpp \
    src/NURBSC.cpp \
    src/NURBSS.cpp \
    src/TRMS.cpp \
    src/BODY.cpp \
    src/DXF_Parser.cpp \
    src/Describe_BODY.cpp \
    src/IGES_Parser.cpp \
    src/KodListFunc.cpp \
    src/MESH.cpp \
    src/Quaternion.cpp \
    src/SFQuant.cpp \
    src/STL_Parser.cpp \
    src/VRML_Parser.cpp \
    src/KodatunoKernel.cpp

HEADERS += \
    include/NURBS.h \
    include/NURBSC.h \
    include/NURBSS.h \
    include/TRMS.h \
    include/BODY.h \
    include/DXF_Parser.h \
    include/Describe_BODY.h \
    include/IGES_Parser.h \
    include/KodListFunc.h \
    include/MESH.h \
    include/Quaternion.h \
    include/SFQuant.h \
    include/STL_Parser.h \
    include/VRML_Parser.h \
    include/KodatunoKernel.h

# Default rules for deployment.
unix {
    target.path = $$[QT_INSTALL_PLUGINS]/generic
}
!isEmpty(target.path): INSTALLS += target
