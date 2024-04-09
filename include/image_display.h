#ifndef IMAGE_DISPLAY_H_GUARD 
#define IMAGE_DISPLAY_H_GUARD 
#include <inttypes.h>
#include <stdbool.h>
typedef void ImageDisplay;
ImageDisplay* create_image_display(const char* display_name,uint32_t display_width,uint32_t display_height);
bool should_close_image_display(ImageDisplay* display);
void refresh_image_display(ImageDisplay* display);
void imshow_image_display(ImageDisplay* display,uint32_t img_width,uint32_t img_height,uint8_t* img_data);
void destroy_image_display(ImageDisplay* display);
#endif
