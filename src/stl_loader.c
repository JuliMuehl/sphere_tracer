#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"


char* read_file(const char* fname,size_t* length){
  FILE* fp = fopen(fname,"r");
  fseek(fp,0,SEEK_END);
  *length= ftell(fp);
  rewind(fp);
  char* text = malloc(sizeof(char) * *length);
  fread(text,1,*length,fp);
  return text;
}

struct TriangleArray load_stl(const char* fname){
  size_t text_len;
  char* text = read_file(fname,&text_len);
  size_t triangles_cap = 256;
  size_t triangles_len = 0;
  struct Triangle* triangles = malloc(sizeof(struct Triangle) * triangles_cap);
  size_t i = 0;
  for(;text[i] != '\n' && i < text_len;i++); // skip solid decleration
  i++; // skip \n
  for(;i<text_len;i++){
    if(text[i] == 'e') break; // reached endfacet statement
    struct Triangle* triangle = &triangles[triangles_len];
    sscanf(&text[i],"facet normal %f %f %f",&triangle->normal.x,&triangle->normal.y,&triangle->normal.z);
    for(;text[i] != '\n' && i < text_len;i++); // skip to end of facet normal line
    i++; // skip \n
    for(;text[i] != '\n' && i < text_len;i++); // skip outer loop line
    i++; // skip \n
    sscanf(&text[i],"vertex %f %f %f",&triangle->p1.x,&triangle->p1.y,&triangle->p1.z);
    for(;text[i] != '\n' && i < text_len;i++); // skip to end of vertex ... line
    i++; // skip \n
    sscanf(&text[i],"vertex %f %f %f",&triangle->p2.x,&triangle->p2.y,&triangle->p2.z);
    for(;text[i] != '\n' && i < text_len;i++); // skip to end of vertex ... line
    i++; // skip \n
    sscanf(&text[i],"vertex %f %f %f",&triangle->p3.x,&triangle->p3.y,&triangle->p3.z);
    for(;text[i] != '\n' && i < text_len;i++); // skip to end of vertex ... line
    i++; // skip \n
    for(;text[i] != '\n' && i < text_len;i++); // skip endloop line
    i++; // skip \n
    for(;text[i] != '\n' && i < text_len;i++); // skip endfacet line
    triangles_len++;
    if(triangles_len == triangles_cap){
      triangles_cap *= 2;
      triangles = reallocarray(triangles,triangles_cap,sizeof(struct Triangle));
    }
    //triangles_len++;
  }
  struct TriangleArray res;
  res.data = triangles;
  res.len = triangles_len;
  return res;
}
