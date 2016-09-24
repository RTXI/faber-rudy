PLUGIN_NAME = faber_rudy_2000

RTXI_INCLUDES=/usr/local/lib/rtxi_includes

HEADERS = faber-rudy.h \
          ${RTXI_INCLUDES}/rtmath.h \
          ${RTXI_INCLUDES}/powfast.hpp \

SOURCES = faber-rudy.cpp \
          ${RTXI_INCLUDES}/rtmath.cpp \
          ${RTXI_INCLUDES}/powfast.cpp \

LIBS =

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
