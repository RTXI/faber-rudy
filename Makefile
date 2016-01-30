PLUGIN_NAME = faber_rudy_2000

HEADERS = faber-rudy.h \
          include/RealTimeMath.h \
          include/PowFast.hpp

SOURCES = faber-rudy.cpp \
          include/RealTimeMath.cpp \
          include/PowFast.cpp

LIBS =

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
