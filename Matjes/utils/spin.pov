#include "colors.inc"
#include "shapes.inc"
#include "metals.inc"
#include "finish.inc"
#include "textures.inc"

#declare DarkGrey = rgb <.4,.4,.4>;

//#declare stretch=1.0;

//global_settings { ambient_light rgb<1., 1., 1.> max_trace_level 20  assumed_gamma 2.2 }

//***************************************
//             Lichtquellen
//***************************************

light_source { <12, -12, -12> 
    color White
shadowless
}
/*light_source {
    <30, 45, 90>
    color White
    spotlight
    radius 15
    falloff 20
    tightness 10
    point_at <10, 10, 20>
  }
*/
light_source { <12, 0, -12> 
    color White
    area_light <5, 0, 0>, <0, 5, 0>, 5, 5
    adaptive 1
    jitter
}
//light_source { <-90, -94, 93> color White}

//***************************************
//               Camera
//***************************************

camera {
  sky   <0,0,1>
//top view//
  location <0,0,200 > //top
  look_at  <50,0,0>
  right    <0,-1,0>
//  location <35, 25 , 34 > // #1
//  look_at  <15, 15  2.>   // #1
//  location <35, 30 ,-70 > 
//  look_at  <35, 36  ,2.> 

  angle 40
}

background { color White }

#declare Spinfont="arial.ttf"

//----------------------------------------------------
//          Reference Frame
//----------------------------------------------------

#macro Spin(theta,phi,px,py,pz,Rc,Bc,Gc)
object{
cone{
     <0.0,0.0,1.0>, 0.0  // Center and radius of one end
     <0.0,0.0,-1.0>, 0.18 // Center and radius of other end
rotate    <0.0,theta,0.0>              // phi rotation
rotate    <0.0,0.0,phi>            // theta rotation
translate <px,py,pz>                 // position of the vector
texture {pigment  {color rgb <Rc,Bc,Gc>}}      //colors

     }}
#end

#macro Vector(length,theta,phi,px,py,pz,Rc,Bc,Gc)
object{
cone{
          <0.0,0.0,length/2.0>, 0.0  // Center and radius of one end
          <0.0,0.0,-length/2.0>, 0.18 // Center and radius of other end
rotate    <0.0,theta,0.0>              // phi rotation
rotate    <0.0,0.0,phi>            // theta rotation
translate <px,py,pz>                 // position of the vector
texture{pigment  {color rgb <Rc,Bc,Gc>}}      //colors

     }}
#end

#macro Density(rho,theta,phi,px,py,pz,Rc,Bc,Gc)
object{
cone{
          <0.0,0.0,1.0>, 0.0  // Center and radius of one end
          <0.0,0.0,-1.0>, 0.18 // Center and radius of other end
rotate    <0.0,phi,0.0>              // phi rotation
rotate    <0.0,0.0,theta>            // theta rotation
translate <px,py,pz>                 // position of the vector
texture{pigment  {color rgb <Rc,Bc,Gc>}}      //colors
     }}
object{
box{
          
}
}
#end

#include "Spinse1.0000.dat"

