CXX=g++
CXXFLAGS=-Wall -O3 -mavx
DEFINES=
INCLUDES=-I../
THIRDLIB=-L../
LIBS=../libsimd.a

TARGET=test_simd
SRCS=$(wildcard *.cpp)
OBJDIR=./objs
OBJS=$(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))
DEPS=$(patsubst %.cpp,$(OBJDIR)/%.d,$(SRCS))

all:$(TARGET)

$(OBJDIR)/%.o:%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MF $(@:.o=.d) -MT $@ -c $< -o $@

-include $(DEPS)

$(TARGET):$(OBJS)
	$(CXX) -o $@ $(OBJS) $(THIRDLIB) $(LIBS)

clean:
	rm -rf $(DEPS) $(OBJS) $(OBJDIR) $(TARGET)
