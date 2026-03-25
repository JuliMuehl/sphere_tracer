CC=gcc
EXECUTABLE=sphere_tracer
LIBS = -lm -lglfw -lGLEW -lGL
OPTIM=-Ofast

OBJS = aabbtree.o\
	   bsdf.o\
	   image_display.o\
	   stl_loader.o\

$(EXECUTABLE): $(OBJS)
	$(CC) $(OPTIM) -Iext -Iinclude -pthread  $(LIBS) $(OBJS) src/sphere_tracer.c -o $(EXECUTABLE) 

image_display.o:
	$(CC) $(OPTIM) -Iinclude -c -pthread  $(LIBS) src/image_display.c -o image_display.o

bsdf.o:
	$(CC) $(OPTIM) -Iinclude -c -pthread  $(LIBS) src/bsdf.c -o bsdf.o

aabbtree.o:
	$(CC) $(OPTIM) -Iinclude -c -pthread  $(LIBS) src/aabbtree.c -o aabbtree.o

stl_loader.o:
	$(CC) $(OPTIM) -Iinclude -c -pthread $(LIBS) src/stl_loader.c -o stl_loader.o

PHONY: clean

clean:
	rm *.o
