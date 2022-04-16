//#############################################//
//######## This is an example script ##########//
//######## For more information, see ##########//
// https://www.povray.org/documentation/3.7.0/ //
//#############################################//

#include "colors.inc"
#include "shapes.inc"
#include "metals.inc"
#include "finish.inc"
#include "textures.inc"
#include "functions.inc"

//##################//
// Background color //
//##################//

background { color rgb < 0.03, 0.03, 0.03 > }	// Defines the background color in < r, g, b >, Here: darkgrey

//################//
// Setup of Light //
//################//

// Example for pointlight //

light_source {
       	<0, 0, 25> 	// Location of the source in < x, y, z >
    	color White 	// Color of the light
}

// Example for Spotlight //

/*light_source { 		// Light Source FLAG
	<0, 0, 12>		// Location of the source in < x, y, z >
	point_at <0, 0, 0>	// Location of the point where light is pointing at in < x, y, z >
    	color rgb < 10.0, 10.0, 10.0 >  // Color of the light in < r, g, b >
	spotlight	// Light Type (spotlight, shadowless, cylinder, parallel)
	radius 35 degrees	// Radius angle of the Spotlight: angle from the center line
	falloff 45 degrees	// Falloff angle of the Spotlight: angle from the center line
	tightness 0
}
*/

//#############################//
// Photon Mapping in the Scene //
//#############################//

/*
photons {
  target
  reflection on
  refraction on
}
*/

//#####################//
// Setup of the Camera //
//#####################//

// Example for easy camera //

camera {


//location < 0, 0, 0 > 	// Location of the camera in < x, y, z >
//look_at < 0, 0, 0 > 	// Position where the camera should look at in < x, y, z >
	//sky < 0.0, 0, 0 > 	// Tilt of the camera to motify the look_at FLAG - very counter-intuitive
	angle 0	// Viewing angel in degrees
	right x
}

/*
Sketch of the applied coordinate system for the cones
        	y
        	^
        	|
        	|
        	|
	x<-------	
*/


//#################//
// Defining macros //
//#################//

#macro Spin(pos_x,pos_y,pos_z,m_x,m_y,m_z)	// definition of the Macro Spin, dependent on 3 coordinates and 3 coordinates for magnetic moment
object{
#if(m_x != 0 | m_y != 0)
	#declare theta = atan2(m_y,m_x)/(2*pi) *360;
#else
	#declare theta = 0;
#end
cone{				// The spins will be visualized as cones
	< 0.0, 0.0, -0.5 >, 0.3	// < Pos of Base point >, Radius Base
	< 0.0, 0.0, +0.5 >, 0.0 // < Pos of Cap point >, Radius Cap
rotate < 0, -(acos((m_z) / sqrt(m_x*m_x + m_y*m_y + m_z*m_z )))/(2*pi) *360 , -theta > // rotation around < x, y, z > axis
translate < -pos_x, pos_y, pos_z > // Position of the vector

/*
Sketch of the applied coordinate system for the cones
	y
	^
	|
	|
	|
	----------> x
*/


// Texture (surface of the object) with colorcode in < r, g, b > //

texture{pigment  {color rgb < (1+m_z) , 1-abs(m_z) , (abs(m_z -1)) >} // Colorcode for the cones: z = strongred, -z = strongblue, 0z = lightgrey/white
}
}

}
#end


#include "position_magnetization.dat"
