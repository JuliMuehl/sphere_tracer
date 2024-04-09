CC=cc
LD=ld
clean:
	rm *.o
sphere_tracer:*.o
	$(CC)  -lm  -pthread -lglfw -lGLEW -lGL -o sphere_tracer *.o
*.o: src/*.c
	$(CC) -Iext -Iinclude -Ofast -c src/*.c 
run:sphere_tracer
	./sphere_tracer 
