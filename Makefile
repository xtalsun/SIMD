CXX=g++
CXXFLAGS=-Wall -O3 -mavx
AR=ar

DEFINES=
INCLUDES=

TARGET=libsimd.a
SRCS=$(wildcard *.cpp)
OBJDIR=./objs
OBJS=$(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))
DEPS=$(patsubst %.cpp,$(OBJDIR)/%.d,$(SRCS))

all:$(TARGET)

$(OBJDIR)/%.o:%.cpp
	#$(@D) == $(OBJDIR)/%.o -> $(OBJDIR)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -MMD -MF $(@:.o=.d) -MT $@ -c -fPIC $< -o $@

-include $(DEPS)

$(TARGET):$(OBJS)
	$(AR) rcs $@ $(OBJS)

clean:
	rm -rf $(DEPS) $(OBJS) $(OBJDIR) $(TARGET)
