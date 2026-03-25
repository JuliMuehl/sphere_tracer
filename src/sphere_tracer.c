#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <float.h>
#include <complex.h>
#include <pthread.h>
#include <unistd.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "vectors.h"
#include "image_display.h"
#include "bsdf.h"
#include "lcg_rng.h"
#include "geom.h"
#include "sampling.h"
#include "mesh.h"
#include "aabbtree.h"
#include "world.h"

struct Vector3f irradiance_bsdf_pathtracer(struct Ray* ray,struct World* world,struct LCGState* rng){
  struct Vector3f f = {0.0f,0.0f,0.0f};
  struct Vector3f throughput = {1.0f,1.0f,1.0f};
  float alpha = 0.99;
  uint32_t i = 0;
  float mean_free_path = 15.0f;
  for(;;){
    struct Intersection it = intersect_ray_world(ray,world);
    if(it.hit && ((dot_v3f(it.normal,ray->direction) >= 0.0f) && it.bsdf.type != BSDF_DIELECTRIC)){
      break;
    }
    bool rr_killed = false;
    while(it.bsdf.type == BSDF_DIELECTRIC){
      if(next_float_lcg(rng) < (1.0-alpha)){
        rr_killed = true;
        break;
      }
      struct Frame frame = frame_from_normal(it.normal);
      struct BSDFQueryRecord rec;
      struct Vector2f sample = {next_float_lcg(rng),next_float_lcg(rng)};
      rec.wi = frame_to_local(&frame,smul_v3f(-1.0f,ray->direction)); sample_bsdf(&it.bsdf,&rec,sample);
      throughput.x /= alpha;
      throughput.y /= alpha;
      throughput.z /= alpha;
      ray->direction = frame_to_global(&frame,rec.wo);
      ray->origin = add_v3f(it.point,smul_v3f(1e-3,ray->direction));
      ray->tmin = 1e-2;
      it = intersect_ray_world(ray,world);
    }
    if(rr_killed) break;
    if(next_float_lcg(rng) < (1.0-alpha)){
      rr_killed = true;
      break;
    }
    throughput.x /= alpha;
    throughput.y /= alpha;
    throughput.z /= alpha;
    float scatter_t = -logf(next_float_lcg(rng)) * mean_free_path;
    struct Vector3f scatter_point = add_v3f(ray->origin,smul_v3f(scatter_t,ray->direction));
    ray->origin = scatter_point;
    if(scatter_t <= it.t){
      struct Vector2f sample = {next_float_lcg(rng),next_float_lcg(rng)};
      float p_sample_uniform_volume = 0.5f;
      bool should_sample_uniform = next_float_lcg(rng) < p_sample_uniform_volume;
      struct Vector3f y = {0.0f,0.0f,0.0f},dist = {1.0f,1.0f,1.0f};
      bool hit_light_source = false;
      if(should_sample_uniform){
        ray->direction = warp_square_to_uniform_hemisphere(sample);
        if(next_float_lcg(rng) > 0.5){
          ray->direction = smul_v3f(-1.0f,ray->direction);
        }  
        struct Intersection it_light = intersect_ray_area_light(ray,&world->light,false);
        hit_light_source = it_light.hit;
        y = it_light.point;
        dist = sub_v3f(y,ray->origin);
      }else{
        y = warp_square_to_area_light(&world->light,sample);
        ray->direction = normalize_v3f(sub_v3f(y,ray->origin));
        hit_light_source = true;
        dist = sub_v3f(y,ray->origin);
      }
      float A = (world->light.xmax - world->light.xmin) * (world->light.zmax - world->light.zmin);
      float dist_len = norm_v3f(dist);
      float J = fabs(dist.y) / (dist_len*dist_len*dist_len);
      if(!hit_light_source) J = 1.0f;
      float pdf = p_sample_uniform_volume + hit_light_source * (1.0f-p_sample_uniform_volume) / J / A * 4.0f * M_PI;
      throughput.x /= pdf;
      throughput.y /= pdf;
      throughput.z /= pdf;
      continue;
    }

    if(it.hit){
      f.x += it.emission.x * throughput.x;
      f.y += it.emission.y * throughput.y;
      f.z += it.emission.z * throughput.z;
      if(it.bsdf.type == BSDF_ZERO){
        break;
      }
      if(next_float_lcg(rng) < (1.0-alpha)){
        break;
      }
      
      float p_sample_bsdf = 0.5;
      struct Frame frame = frame_from_normal(it.normal);
      struct BSDFQueryRecord rec;
      struct Vector2f sample = {next_float_lcg(rng),next_float_lcg(rng)};
      struct Vector3f x = it.point,y = warp_square_to_area_light(&world->light,sample);
      rec.wi = frame_to_local(&frame,smul_v3f(-1.0f,ray->direction));
      struct Vector3f dist = sub_v3f(y,x);
      rec.wo = frame_to_local(&frame,normalize_v3f(dist));
      float hit_light_source = 1.0f;
      bool surface_is_glossy = it.bsdf.type == BSDF_MIRROR ||
        (it.bsdf.type == BSDF_COOK_TORRANCE && it.bsdf.cook_torrance.alpha < 0.1f);
      if(surface_is_glossy){
        p_sample_bsdf = 1.0f;
      }
      bool should_sample_bsdf = next_float_lcg(rng) < p_sample_bsdf;
      if(should_sample_bsdf){
        sample_bsdf(&it.bsdf,&rec,sample);
        struct Ray ray_light = {1e-5f,x,frame_to_global(&frame,rec.wo)};
        struct Intersection it_light = intersect_ray_area_light(&ray_light,&world->light,false);
        hit_light_source = (float)(it_light.hit);
        y = it_light.point;
      }
      struct Vector3f bsdf = eval_bsdf(&it.bsdf,&rec);
      float pdf = p_sample_bsdf * pdf_bsdf(&it.bsdf,&rec);
      if(should_sample_bsdf && (it.bsdf.type == BSDF_MIRROR || it.bsdf.type == BSDF_DIELECTRIC)){
        pdf = p_sample_bsdf * rec.wo.z;
      }
      float A = (world->light.xmax - world->light.xmin) * (world->light.zmax - world->light.zmin);
      float dist_len = norm_v3f(dist);
      float J = fabs(dist.y) / (dist_len*dist_len*dist_len);
      
      if(hit_light_source == 0.0f) J = 1.0f;
      pdf += (1.0f - p_sample_bsdf) * hit_light_source / A / J;
      if(isnan(pdf) || isnan(1.0f/pdf)){
        throughput = smul_v3f(0.0f,throughput);
        break;
      }
      if(rec.wi.z > 0.0f && rec.wo.z > 0.0f){
        throughput.x *= bsdf.x * rec.wo.z / pdf / alpha;
        throughput.y *= bsdf.y * rec.wo.z / pdf / alpha;
        throughput.z *= bsdf.z * rec.wo.z / pdf / alpha;
      }else{
        throughput.x *= 0.0f;
        throughput.y *= 0.0f;
        throughput.z *= 0.0f;
      }
      
      ray->direction = frame_to_global(&frame,rec.wo);
      ray->origin = add_v3f(it.point,smul_v3f(1e-3,ray->direction));
      ray->tmin = 1e-2;
    }else{
      f.x += world->sky_color.x*throughput.x;
      f.y += world->sky_color.y*throughput.y;
      f.z += world->sky_color.z*throughput.z;
      break;
    }
  }
  return f;
}

struct Vector3f irradiance_next_event_estimation_area_light_only(struct Ray* ray,struct World* world,struct LCGState* rng){
  struct Vector3f f = {0.0f,0.0f,0.0f};
  struct Vector3f throughput = {1.0f,1.0f,1.0f};
  float alpha = 0.99;
  struct Intersection it = intersect_ray_world(ray,world);
  struct Intersection it_init = it;
  if(it.hit){
    if(it.bsdf.type == BSDF_ZERO){
      return it.emission;
    }
  }
  for(;;){
    if(it.hit && dot_v3f(it.normal,ray->direction) >= 0.0f){
      break;
    }
    bool rr_killed = false;
    bool from_dielectric = false;
    while(it.bsdf.type == BSDF_DIELECTRIC){
      from_dielectric = true;
      if(next_float_lcg(rng) < (1.0-alpha)){
        rr_killed = true;
        break;
      }
      struct Frame frame = frame_from_normal(it.normal);
      struct BSDFQueryRecord rec;
      struct Vector2f sample = {next_float_lcg(rng),next_float_lcg(rng)};
      rec.wi = frame_to_local(&frame,smul_v3f(-1.0f,ray->direction));
      sample_bsdf(&it.bsdf,&rec,sample);
      throughput.x /= alpha;
      throughput.y /= alpha;
      throughput.z /= alpha;
      ray->direction = frame_to_global(&frame,rec.wo);
      ray->origin = it.point;
      ray->tmin = 1e-2;
      it = intersect_ray_world(ray,world);
    }
    if(rr_killed) break;
    if(from_dielectric && it.hit && it.bsdf.type == BSDF_ZERO){
      f.x += throughput.x  * it.emission.x;
      f.y += throughput.y  * it.emission.y;
      f.z += throughput.z  * it.emission.z;
      break;
    }
    if(it.hit){
      if(it.bsdf.type == BSDF_ZERO){
        break;
      }
    
      struct Frame frame = frame_from_normal(it.normal);
      struct BSDFQueryRecord rec;
      struct Vector2f sample = {next_float_lcg(rng),next_float_lcg(rng)};
      rec.wi = frame_to_local(&frame,smul_v3f(-1.0f,ray->direction));
      //Next event estimation
      struct Vector3f x = it.point,y = warp_square_to_area_light(&world->light,sample);
      struct Vector3f dist = sub_v3f(y,x);
      float dist_len = norm_v3f(dist);
      struct Ray ray_light;
      ray_light.tmin = 1e-2f;
      ray_light.origin = x;
      ray_light.direction = normalize_v3f(dist);
      struct Intersection it_light = intersect_ray_world(&ray_light,world);
      float J = fabs(dist.y) / (dist_len*dist_len*dist_len);
      float A = (world->light.xmax - world->light.xmin) * (world->light.zmax - world->light.zmin);
      if(it_light.hit){
        struct BSDFQueryRecord rec_direct_light;
        rec_direct_light.wi = rec.wi;
        rec_direct_light.wo = smul_v3f(1.0f,frame_to_local(&frame,ray_light.direction));
        struct Vector3f bsdf = eval_bsdf(&it.bsdf,&rec_direct_light);
        float pdf = 1.0f / A / J;
        f.x += throughput.x * bsdf.x * fabs(rec_direct_light.wo.z) / pdf  * it_light.emission.x;
        f.y += throughput.y * bsdf.y * fabs(rec_direct_light.wo.z) / pdf  * it_light.emission.y;
        f.z += throughput.z * bsdf.z * fabs(rec_direct_light.wo.z) / pdf  * it_light.emission.z;
      }
      if(next_float_lcg(rng) < (1.0-alpha)){
        break;
      }
      //BSDF sampling
      sample = (struct Vector2f){next_float_lcg(rng),next_float_lcg(rng)};
      sample_bsdf(&it.bsdf,&rec,sample);
      struct Vector3f bsdf = eval_bsdf(&it.bsdf,&rec);
      float pdf = pdf_bsdf(&it.bsdf,&rec);

      throughput.x *= bsdf.x * fabs(rec.wo.z) / pdf / alpha;
      throughput.y *= bsdf.y * fabs(rec.wo.z) / pdf / alpha;
      throughput.z *= bsdf.z * fabs(rec.wo.z) / pdf / alpha;
      ray->direction = frame_to_global(&frame,rec.wo);
      ray->origin = it.point;
      ray->tmin = 1e-2;
      it = intersect_ray_world(ray,world);
    }else{
      f.x += world->sky_color.x*throughput.x;
      f.y += world->sky_color.y*throughput.y;
      f.z += world->sky_color.z*throughput.z;
      break;
    }
  }
  return f;
}

struct WorkerData{
  uint32_t imin,jmin,imax,jmax;
  struct World* world;
  uint32_t num_rays;
  uint32_t img_width,img_height;
  uint8_t* img_data;
  struct KDTreeNode* photon_map;
};

void* render_job(void* data){
  struct WorkerData* worker_data = data;
  struct World* world = worker_data->world;
  uint32_t imin = worker_data->imin,jmin = worker_data->jmin,imax = worker_data->imax,jmax = worker_data->jmax;
  uint8_t* img_data = worker_data->img_data;
  uint32_t img_width = worker_data->img_width,img_height = worker_data->img_height;
  uint32_t num_rays = worker_data->num_rays;
  struct LCGState rng;
  seed_lcg(&rng,(((((((imin+5381) << 5) + imax) << 5) + jmin) << 5) + jmax ) << 5);
  float aspect = (float) img_width / img_height;
  for(uint32_t i = imin;i<imax;i++){
    for(uint32_t j = jmin;j<jmax;j++){
      float u=(float)j/img_width,v=1.0f-(float)i/img_height;
      struct Ray init_ray;
      init_ray.direction = normalize_v3f((struct Vector3f){0.01 + 2.0f * u - 1.0f,2.0f * v - 1.0f,1.0f});
      init_ray.direction.y /= aspect;
      init_ray.origin = (struct Vector3f) {-0.0,1.0f,0.5f};
      init_ray.tmin = 0.0f;
      struct Intersection it = intersect_ray_world(&init_ray,world);
      struct Vector3f color = {0.0f,0.0f,0.0f};
      bool volumetric = true;
      if(!volumetric && !it.hit){
        color = world->sky_color;
      }else if(!volumetric && it.bsdf.type == BSDF_ZERO){
        color = it.emission;
      }else{
        for(uint32_t ray_idx = 0;ray_idx<num_rays;ray_idx++){
          struct Ray ray = init_ray;
          //color = add_v3f(color,irradiance_bsdf_pathtracer(&ray,world,&rng));
          color = add_v3f(color,irradiance_next_event_estimation_area_light_only(&ray,world,&rng));
        }
        color = smul_v3f(1.0f/num_rays,color);
      }
      color.x = fmax(fmin(color.x,1.0f),0.0f);
      color.y = fmax(fmin(color.y,1.0f),0.0f);
      color.z = fmax(fmin(color.z,1.0f),0.0f);
      img_data[i*img_width*3 + j*3 + 0] = (uint8_t)(color.x * 255);
      img_data[i*img_width*3 + j*3 + 1] = (uint8_t)(color.y * 255);
      img_data[i*img_width*3 + j*3 + 2] = (uint8_t)(color.z * 255);
    }
  }
  return NULL;
}

struct DisplayThreadData{
  uint32_t img_width,img_height;
  uint8_t* img_data;
};

void* display_job(void* data){
  struct DisplayThreadData* thread_data = (struct DisplayThreadData*) data;
  uint32_t img_width=thread_data->img_width,img_height = thread_data->img_height;
  uint8_t* img_data = thread_data->img_data;
  ImageDisplay* display = create_image_display(__FILE__,img_width,img_height);
  while(!should_close_image_display(display)){
    imshow_image_display(display,img_width,img_height,img_data);
    refresh_image_display(display);
  }
  destroy_image_display(display);
  exit(0);
  return NULL;
}

void render_multithreaded(const char* img_path,uint32_t img_width,uint32_t img_height,struct World* world,uint32_t num_rays,uint32_t num_threads){
  uint8_t* img_data = malloc(sizeof(uint8_t) * img_width * img_height * 3);
  memset(img_data,0,img_width*img_height*3*sizeof(uint8_t));
  struct WorkerData* worker_data = malloc(sizeof(worker_data[0])*num_threads);
  pthread_t* threads = malloc(sizeof(threads[0])*num_threads);
  struct DisplayThreadData display_thread_data;
  display_thread_data.img_width = img_width;
  display_thread_data.img_height = img_height;
  display_thread_data.img_data = img_data;
  pthread_t display_thread;
  pthread_create(&display_thread,NULL,display_job,&display_thread_data);
  
  for(uint32_t i = 0;i<num_threads;i++){
    worker_data[i].img_width = img_width;
    worker_data[i].img_height = img_height;
    worker_data[i].img_data = img_data;
    worker_data[i].num_rays = num_rays;
    worker_data[i].world = world;
    worker_data[i].imin = (img_height/num_threads) * i;
    worker_data[i].imax = (img_height/num_threads) * (i+1);
    worker_data[i].jmin = 0;
    worker_data[i].jmax = img_width;
  }

  worker_data[num_threads-1].imax = img_height;
  for(uint32_t i = 0;i<num_threads;i++){
    int err = pthread_create(&threads[i],NULL,render_job,&worker_data[i]);
    if(err) printf("Error running pthread_create for threads[%d]\n",i);
  }

  for(uint32_t i = 0;i<num_threads;i++){
    int err = pthread_join(threads[i],NULL);
    if(err) printf("Error running pthread_join for threads[%d]\n",i);
  }
  
  stbi_write_png(img_path,img_width,img_height,3,&img_data[0],0);
  free(threads);
  free(worker_data);
  int err = pthread_join(display_thread,NULL);
  free(img_data);
}

void render(const char* img_path,uint32_t img_width,uint32_t img_height,struct World* world,uint32_t num_rays){
  uint8_t* img_data = malloc(sizeof(uint8_t) * img_width * img_height * 3);
  struct LCGState rng;
  seed_lcg(&rng,1);
  float aspect = (float) img_width / img_height;
  for(int i = 0;i<img_height;i++){
    for(int j = 0;j<img_width;j++){
      float u=(float)j/img_width,v=1.0f-(float)i/img_height;
      struct Ray init_ray;
      init_ray.direction = normalize_v3f((struct Vector3f){2.0f * u - 1.0f,2.0f * v - 1.0f,1.0f});
      init_ray.direction.y /= aspect;
      init_ray.origin = (struct Vector3f) {-0.0,1.0f,0.5f};
      init_ray.tmin = 0.0f;

      struct Vector3f color = {0.0f,0.0f,0.0f}; // constant environment light with intensity 1.0f
      for(uint32_t ray_idx = 0;ray_idx<num_rays;ray_idx++){
        struct Ray ray = init_ray;
        color = add_v3f(color,irradiance_bsdf_pathtracer(&ray,world,&rng));
        //color = add_v3f(color,irradiance_next_event_estimation_area_light_only(&ray,world,&rng));
      }
      color = smul_v3f(1.0f/num_rays,color);
      
      color.x = fmax(fmin(color.x,1.0f),0.0f);
      color.y = fmax(fmin(color.y,1.0f),0.0f);
      color.z = fmax(fmin(color.z,1.0f),0.0f);
      img_data[i*img_width*3 + j*3 + 0] = (uint8_t)(color.x * 255);
      img_data[i*img_width*3 + j*3 + 1] = (uint8_t)(color.y * 255);
      img_data[i*img_width*3 + j*3 + 2] = (uint8_t)(color.z * 255);
    }
  }

  stbi_write_png(img_path,img_width,img_height,3,&img_data[0],0);
  free(img_data);
}

int main(int argc,char** argv){
  struct TriangleArray triangles = load_stl("/home/juli/sphere_tracer/monkey.stl");
  struct World world;
  world.aabbtree = create_aabbtree(triangles);
  set_aabbtree_position(world.aabbtree,(struct Vector3f){0.0f,1.5f,7.0f});
  float light_x = 0.0;
  float light_z = 2.0;
  world.light.y = 3.5;
  world.light.xmin = -1.5 + light_x;
  world.light.xmax =  1.5 + light_x;
  world.light.zmin =  -0.5 + light_z;
  world.light.zmax =  0.5 + light_z;
  world.sky_color = (struct Vector3f){0.02f,0.02f,0.02f};
  float intensity = 25.0;
  struct Vector3f light_color = (struct Vector3f){1.0,1.0,1.0};
  world.light.emission_intensity = smul_v3f(intensity,light_color);
  world.num_spheres = 6;
  world.spheres = malloc(sizeof(world.spheres[0]) * world.num_spheres);
  world.spheres[0].radius = 1.0f;
  world.spheres[0].center = (struct Vector3f) {-3.0f,1.2f,5.0f};

  world.spheres[1].radius = 1.0f;
  world.spheres[1].center = (struct Vector3f) { 0.0f,1.2f,3.0f};

  world.spheres[2].radius = 1.0f;
  world.spheres[2].center = (struct Vector3f) { 3.0f,1.2f,5.0f};

  world.spheres[3].radius = 1.0f;
  world.spheres[3].center = (struct Vector3f) { 0.0f,1000.2f,8.0f};

  world.spheres[4].radius = 1.0f;
  world.spheres[4].center = (struct Vector3f) { 3.0f,1.2f,8.0f};
  
  world.spheres[5].radius = 1.0f;
  world.spheres[5].center = (struct Vector3f) {-3.0f,1.2f,8.0f};

  for(uint32_t i = 0;i<world.num_spheres;i++){
    world.spheres[i].emission = (struct Vector3f) {0.0f,0.0f,0.0f};
  }

  world.plane_y = 0.0f;

  struct BSDF red_bsdf;
  struct BSDF blue_bsdf;
  struct BSDF mirror_bsdf;
  red_bsdf.type = BSDF_DIFFUSE;
  red_bsdf.diffuse.albedo = (struct Vector3f){0.5,0.2,0.2};
  blue_bsdf.type = BSDF_DIFFUSE;
  blue_bsdf.diffuse.albedo = (struct Vector3f){0.2,0.3,0.4};
  mirror_bsdf.type = BSDF_MIRROR;

  struct BSDF copper_bsdf;
  struct BSDF silver_bsdf;
  struct BSDF gold_bsdf;
  gold_bsdf.type = BSDF_COOK_TORRANCE;
  gold_bsdf.cook_torrance.eta = (struct Vector3f){0.1715, 0.382257, 1.43639};
  gold_bsdf.cook_torrance.kappa = (struct Vector3f){3.97557, 2.38041, 1.6045};

  silver_bsdf.type = BSDF_COOK_TORRANCE;
  silver_bsdf.cook_torrance.eta = (struct Vector3f){0.154804, 0.116449, 0.138004};
  silver_bsdf.cook_torrance.kappa= (struct Vector3f){4.82855, 3.12125, 2.14529};

  copper_bsdf.type = BSDF_COOK_TORRANCE;
  copper_bsdf.cook_torrance.eta = (struct Vector3f){0.211669, 0.919456, 1.10127};
  copper_bsdf.cook_torrance.kappa= (struct Vector3f){3.92127, 2.45671, 2.13517};
  struct BSDF dielectric_bsdf;
  dielectric_bsdf.type = BSDF_DIELECTRIC;
  dielectric_bsdf.dielectric.eta_inside = 1.8;
  dielectric_bsdf.dielectric.eta_outside = 1.0;
  float alpha = 0.10;
  copper_bsdf.cook_torrance.alpha = alpha;
  silver_bsdf.cook_torrance.alpha = alpha;
  gold_bsdf.cook_torrance.alpha = alpha;
  world.spheres[0].bsdf = silver_bsdf;
  world.spheres[1].bsdf = dielectric_bsdf;
  world.spheres[2].bsdf = copper_bsdf; 
  world.spheres[3].bsdf = blue_bsdf;
  world.spheres[4].bsdf = red_bsdf;
  world.spheres[5].bsdf = red_bsdf;
  world.aabbtree_bsdf = gold_bsdf;
  world.plane_bsdf.diffuse.albedo = (struct Vector3f){0.2,0.2,0.2};
  world.plane_bsdf.type = BSDF_DIFFUSE;
  int num_rays = argc >= 2 ? atoi(argv[1]):40;
  //render_multithreaded("out.png",512,288,&world,num_rays,25);
  render_multithreaded("out.png",640,480,&world,num_rays,25);
  //render("out.png",512,288,&world,num_rays);
  free(world.spheres);
}
