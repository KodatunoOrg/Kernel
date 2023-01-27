CONFIG -= qt

TEMPLATE = lib
CONFIG += staticlib

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += include/

SOURCES += \
    src/BODY.cpp \
    src/DXF_Parser.cpp \
    src/Describe_BODY.cpp \
    src/IGES_Parser.cpp \
    src/KodListFunc.cpp \
    src/MESH.cpp \
    src/NURBS_Func.cpp \
    src/Quaternion.cpp \
    src/SFQuant.cpp \
    src/STL_Parser.cpp \
    src/StdAfx.cpp \
    src/VRML_Parser.cpp

HEADERS += \
    include/BODY.h \
    include/DXF_Parser.h \
    include/Describe_BODY.h \
    include/IGES_Parser.h \
    include/KodListFunc.h \
    include/MESH.h \
    include/NURBS_Func.h \
    include/Quaternion.h \
    include/SFQuant.h \
    include/STL_Parser.h \
    include/StdAfx.h \
    include/VRML_Parser.h

# Default rules for deployment.
unix {
    target.path = $$[QT_INSTALL_PLUGINS]/generic
}
!isEmpty(target.path): INSTALLS += target
