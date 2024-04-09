#include "image_display.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#define STRING(x) #x

GLint compile_shader(GLenum type,const char* source){
  GLint shader = glCreateShader(type);
  GLint source_len = strlen(source);
  glShaderSource(shader,1,&source,&source_len);
  glCompileShader(shader);
  GLint status;
  glGetShaderiv(shader,GL_COMPILE_STATUS,&status);
  if(status != GL_TRUE){
    GLint info_log_length;
    glGetShaderiv(shader,GL_INFO_LOG_LENGTH,&info_log_length);
    char* info_log = malloc(sizeof(info_log[0]) * info_log_length);
    glGetShaderInfoLog(shader,info_log_length,NULL,info_log);
    printf("%s shader compile error: %s\n",type == GL_VERTEX_SHADER ? "vertex":"fragment",info_log);
    free(info_log);
    exit(1);
  }
  return shader;
}

GLint compile_program(GLint vertex_shader,GLint fragment_shader){
  GLint program = glCreateProgram();
  glAttachShader(program,vertex_shader);
  glAttachShader(program,fragment_shader);
  glLinkProgram(program);
  GLint status;
  glGetProgramiv(program,GL_LINK_STATUS,&status);
  if(status != GL_TRUE){
    GLint info_log_length;
    glGetProgramiv(program,GL_INFO_LOG_LENGTH,&info_log_length);
    char* info_log = malloc(sizeof(info_log[0]) * info_log_length);
    glGetProgramInfoLog(program,info_log_length,NULL,info_log);
    printf("program link error: %s\n",info_log);
    free(info_log);
    exit(1);
  }
  return program;
}

GLuint load_texture(GLuint texture,int width,int height,uint8_t* img_data){
  glBindTexture(GL_TEXTURE_2D,texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,width,height,0,GL_RGB,GL_UNSIGNED_BYTE,img_data);
  return texture;
}

struct ImageDisplay{
  GLFWwindow* window;
  GLuint uv_buffer,texture;
  GLint program;
  GLint in_uv_location;
};

ImageDisplay* create_image_display(const char* display_name,uint32_t display_width,uint32_t display_height){
  if(!glfwInit()){
    printf("Error initializing GLFW\n");
    exit(1);
  }
  glfwWindowHint(GLFW_RESIZABLE,GLFW_FALSE);
  GLFWwindow* window = glfwCreateWindow(display_width,display_height,display_name,NULL,NULL);
  glfwMakeContextCurrent(window);
  if(glewInit()){
    printf("Error initializing GLEW\n");
    exit(1);
  }
  
  float uv_data[] = {0.0f,0.0f,0.0f,1.0f,1.0f,1.0f,0.0f,0.0f,1.0f,0.0f,1.0f,1.0f};
  GLuint uv_buffer;
  glGenBuffers(1,&uv_buffer);
  glBindBuffer(GL_ARRAY_BUFFER,uv_buffer);
  glBufferData(GL_ARRAY_BUFFER,sizeof(uv_data),&uv_data[0],GL_STATIC_DRAW);

  const char vertex_shader_source[] = "#version 330\n" STRING(
    in vec2 in_uv;
    out vec2 frag_uv;
    void main(){
      frag_uv = in_uv;
      gl_Position = vec4(2.0 * in_uv - 1.0,0.0,1.0);
    }
  );

  const char fragment_shader_source[] = "#version 330\n" STRING(
    in vec2 frag_uv;
    out vec4 frag_color;
    uniform sampler2D img;
    void main(){
      frag_color = texture(img,vec2(frag_uv.x,1.0-frag_uv.y));
    }
  );

  GLint vertex_shader = compile_shader(GL_VERTEX_SHADER,&vertex_shader_source[0]);
  GLint fragment_shader = compile_shader(GL_FRAGMENT_SHADER,&fragment_shader_source[0]);
  GLint program = compile_program(vertex_shader,fragment_shader);
  glDeleteShader(vertex_shader);
  glDeleteShader(fragment_shader);
  glUseProgram(program);
  GLint in_uv_location = glGetAttribLocation(program,"in_uv");
  glEnableVertexAttribArray(in_uv_location);
  glVertexAttribPointer(in_uv_location,2,GL_FLOAT,GL_FALSE,0,NULL);
  GLuint texture;
  glGenTextures(1,&texture);
  struct ImageDisplay* display = malloc(sizeof(*display));
  display->window = window;
  display->program = program;
  display->in_uv_location = in_uv_location;
  display->uv_buffer = uv_buffer;
  display->texture = texture;
  return display;
}

bool should_close_image_display(ImageDisplay* display_){
  struct ImageDisplay* display = display_;
  return glfwWindowShouldClose(display->window);
}

void refresh_image_display(ImageDisplay* display_){
  struct ImageDisplay* display = display_;
  glfwSwapBuffers(display->window);
  glfwPollEvents();
}

void imshow_image_display(ImageDisplay* display_,uint32_t img_width,uint32_t img_height,uint8_t* img_data){
  struct ImageDisplay* display = display_;
  load_texture(display->texture,img_width,img_height,img_data);
  glUseProgram(display->program);
  glBindBuffer(GL_ARRAY_BUFFER,display->uv_buffer);
  glEnableVertexAttribArray(display->in_uv_location);
  glVertexAttribPointer(display->in_uv_location,2,GL_FLOAT,GL_FALSE,0,NULL);
  glClearColor(0.0f,0.0f,0.0f,1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  glDrawArrays(GL_TRIANGLES,0,6);
}

void destroy_image_display(ImageDisplay* display_){
  struct ImageDisplay* display = display_;
  glDeleteTextures(1,&display->texture);
  glDeleteBuffers(1,&display->uv_buffer);
  glDeleteProgram(display->program);
  glfwDestroyWindow(display->window);
  glfwTerminate();
  free(display);
}
