//------------------------------------------------
//
//  cfd_paint
//
//
//-------------------------------------------------

//-------------------------------------------------
//
//  usage:
//
//  cfd_paint is an interactive paint program
//  in which the user paints density, color, and
//  or divergence sources that flow using
//  computational fluid dynamics and react with 
//  obstructions in the space.
//
//  There are two paint modes.  Typing 'o' puts the
//  program in obstruction painting mode. When you
//  hold down the left mouse button and paint, you
//  will see a black obstruction painted.  This 
//  obstruction may be any shape.
//
//  Typing 's' puts the program in source painting 
//  mode.  Now painting with the left mouse button
//  down injects density into the simulation.
//  The flow it produces evolves as you
//  continue to paint.  The flow bounces off any
//  obstructions that have been painted or are
//  subsequently painted.
//
//  Typing 'b' clears all obstructions, flow, density,
//  and color.
//
//  Typing '=' and '-' brightens and darkens the display.
//
//  Pressing the spacebar starts and stops the flow 
//  evolution. While the evolution is stopped, you
//  can continue painting obstructions.
//
//
//
//-------------------------------------------------




#include <cmath>
#include "CmdLineFind.h"


#ifdef __APPLE__
  #include <OpenGL/gl.h>   // OpenGL itself.
  #include <OpenGL/glu.h>  // GLU support library.
  #include <GLUT/glut.h>
#else
  #include <GL/gl.h>   // OpenGL itself.
  #include <GL/glu.h>  // GLU support library.
  #include <GL/glut.h> // GLUT support library.
#endif

#include <iostream>
#include <OpenImageIO/imageio.h>

#include "cfd.h"

using namespace std;
using namespace lux;
OIIO_NAMESPACE_USING

cfd *cfd_demo;
float *display_map;
float *color_map;
float *obstruction_map;
float *source_map;
float *divergence_source_map;

bool capture_screen;
bool is_on;
bool has_boundary;

int frame;
int iwidth, iheight, size;

float lx, ly;
int nx, ny;
int gs, iop;
int display_mode;
enum 
{
   DISPLAY, VELOCITY, PRESSURE, DENSITY
};

string captured_file_basename;

int paint_mode;
enum{ PAINT_OBSTRUCTION, PAINT_SOURCE, PAINT_DIVERGENCE, PAINT_COLOR };

float scaling_factor;

float h = 0.2; // timestep
float grav = 9.8;

#define BRUSH_SIZE 11
int demo = 1;
float obstruction_brush[BRUSH_SIZE][BRUSH_SIZE];
float source_brush[BRUSH_SIZE][BRUSH_SIZE];
float divergence_brush[BRUSH_SIZE][BRUSH_SIZE];

int xmouse_prev, ymouse_prev;

////////  OpenImageIO reader

void readOIIOImage( const char* fname, float* img  )
{
   int xres, yres, channels;
   ImageInput *in = ImageInput::create (fname);
   if (! in) {return;}
   ImageSpec spec;
   in->open (fname, spec);
   xres = spec.width;
   yres = spec.height;
   channels = spec.nchannels;
   float* pixels = new float[xres*yres*channels];
   in->read_image (TypeDesc::FLOAT, pixels);
   long index = 0;
   for( int j=0;j<yres;j++)
   {
      for( int i=0;i<xres;i++ )
      {
         for( int c=0;c<channels;c++ )
         {
	        img[ (i + xres*(yres - j - 1))*channels + c ] = pixels[index++];
         }
      }
   }

   in->close ();
   delete in;
}

void writeOIIOImage( const char* fname, const float* img )
{
   float* pixels = new float[3*iwidth*iheight];
   long index = 0;
   for( int j=0;j<iheight;j++)
   {
      for( int i=0;i<iwidth;i++ )
      {
         for( int c=0;c<3;c++ )
         {
	        pixels[ (i + iwidth*(iheight - j - 1))*3 + c ] = img[index++];
         }
      }
   }

   ImageOutput *out = ImageOutput::create (fname); 
   if( !out )
   {
      cout << "Not able to write an image to file " << fname << endl;
   }
   else
   {
      ImageSpec spec (iwidth, iheight, 3, TypeDesc::FLOAT); 
      spec.attribute("user", "jtessen");
      spec.attribute("writer", "writeOIIOImage" );
      out->open (fname, spec);
      out->write_image (TypeDesc::FLOAT, pixels);
      out->close (); 
      delete out;
   }
   delete[] pixels;
}

void writeTestImage( const char* fname, const float* img )
{
   float* pixels = new float[iwidth*iheight];
   long index = 0;
   for( int j=0;j<iheight;j++)
   {
      for( int i=0;i<iwidth;i++ )
      {
        
           pixels[ (i + iwidth*(iheight - j - 1)) ] = img[index++];
         
      }
   }

   ImageOutput *out = ImageOutput::create (fname); 
   if( !out )
   {
      cout << "Not able to write an image to file " << fname << endl;
   }
   else
   {
      ImageSpec spec (iwidth, iheight, 1, TypeDesc::FLOAT); 
      spec.attribute("user", "ywang");
      spec.attribute("writer", "writeOIIOImage" );
      out->open (fname, spec);
      out->write_image (TypeDesc::FLOAT, pixels);
      out->close (); 
      delete out;
   }
   delete[] pixels;
}

//--------------------------------------------------------
//
//  Initialization routines
//
//  
// Initialize all of the fields to zero
void Initialize( float *data, int size, float value )
{
#pragma omp parallel for
   for(int i=0;i<size;i++ ) { data[i] = value; }
}

void InitializeBrushes()
{
   int brush_width = (BRUSH_SIZE-1)/2;
   for( int j=-brush_width;j<=brush_width;j++ )
   {
      int jj = j + brush_width;
      float jfactor =  (float(brush_width) - fabs(j) )/float(brush_width);
      for( int i=-brush_width;i<=brush_width;i++ )
      {
         int ii = i + brush_width;
         float ifactor =  (float(brush_width) - fabs(i) )/float(brush_width);
         float radius = (jfactor*jfactor + ifactor*ifactor)/2.0;
         source_brush[ii][jj] = pow(radius,0.5);
         divergence_brush[ii][jj] = pow(radius, 0.5);
         obstruction_brush[ii][jj] = 1.0 - pow(radius, 1.0/4.0);
      }
   }
}

void setNbCores( int nb )
{
  // omp_set_num_threads( nb );
}


//----------------------------------------------------

void ConvertToDisplay(int type)
{ 
    if(has_boundary){
        color_map = cfd_demo->getColorMapWithBoundary();
    }else{
        color_map = cfd_demo->getColorMap();
    }
   
    if(type == DISPLAY){

      for( int j=0;j<iheight;j++ )
        {
#pragma omp parallel for
            for(int i=0;i<iwidth;i++ )
                {
                  int index = i + iwidth*j;
                  float r,g,b;
                  r = color_map[index*3+0];
                  g = color_map[index*3+1];
                  b = color_map[index*3+2];
                  display_map[3*index+0] = r * scaling_factor; 
                  display_map[3*index+1] = g * scaling_factor; 
                  display_map[3*index+2] = b * scaling_factor; 
                }
         } 
   }
 
}


//------------------------------------------
//
//  Painting and display code
//

void resetScaleFactor( float amount )
{
   scaling_factor *= amount;
}

void DabSomePaint( int x, int y )
{
   int brush_width = (BRUSH_SIZE-1)/2;
   int xstart = x - brush_width;
   int ystart = y - brush_width;
   if( xstart < 0 ){ xstart = 0; }
   if( ystart < 0 ){ ystart = 0; }

   int xend = x + brush_width;
   int yend = y + brush_width;
   if( xend >= iwidth ){ xend = iwidth-1; }
   if( yend >= iheight ){ yend = iheight-1; }

   source_map = cfd_demo->getSourceMap();
   obstruction_map = cfd_demo->getObstructionMap();
   divergence_source_map = cfd_demo->getDivergenceSourceMap();

   if( paint_mode == PAINT_OBSTRUCTION )
   {
      for(int ix=xstart;ix <= xend; ix++)
      {
         for( int iy=ystart;iy<=yend; iy++)
	     {
          int index = ix + iwidth*(iheight-iy-1);
          obstruction_map[index] *= obstruction_brush[ix-xstart][iy-ystart];
	    }
      }
   }
   else if( paint_mode == PAINT_SOURCE )
   {
      for(int ix=xstart;ix <= xend; ix++)
      {
         for( int iy=ystart;iy<=yend; iy++)
	     {
            int index = ix + iwidth*(iheight-iy-1);
            source_map[index] += source_brush[ix-xstart][iy-ystart];

	     }
      }
   }
   else if( paint_mode == PAINT_DIVERGENCE)
   {
      for(int ix=xstart;ix <= xend; ix++)
      {
         for( int iy=ystart;iy<=yend; iy++)
       {
            int index = ix + iwidth*(iheight-iy-1);
            divergence_source_map[index] += divergence_brush[ix-xstart][iy-ystart];

       }
      }
   }
   return; 
}

//----------------------------------------------------
//
//  GL and GLUT callbacks
//
//----------------------------------------------------

void cbDisplay( void )
{
   glClear(GL_COLOR_BUFFER_BIT );
   if(display_mode == DISPLAY){
      glDrawPixels( iwidth, iheight, GL_RGB, GL_FLOAT, display_map );
   }
   
   glutSwapBuffers();
   if( capture_screen )
   {
      std::stringstream os; os<<frame;
      string dispframe = os.str();
      if( frame < 1000 ){ dispframe = "0" + dispframe; }
      if( frame < 100 ){ dispframe = "0" + dispframe; }
      if( frame < 10 ){ dispframe = "0" + dispframe; }
      string fname = captured_file_basename + "." + dispframe + ".exr";
      writeOIIOImage( fname.c_str(), display_map ); 
      cout << "Frame written to file " << fname << endl;
   }
}

// animate and display new result
void cbIdle()
{
   if(is_on){
      cfd_demo->updateFluid();
   } 

   ConvertToDisplay(display_mode);
   glutPostRedisplay();
   frame++;
}

void cbOnKeyboard( unsigned char key, int x, int y )
{
   switch (key) 
   {
      case '-': case '_':
      resetScaleFactor( 0.9 );
      break;

      case '+': case '=':
      resetScaleFactor( 1.0/0.9 );
      break;

      case 'r':
      scaling_factor = 1.0;
      break;

      case 'c':
      capture_screen = !capture_screen;

      case 'o':
      paint_mode = PAINT_OBSTRUCTION;
      break;

      case 's':
      paint_mode = PAINT_SOURCE;
      break;

      case 'd':
      paint_mode = PAINT_DIVERGENCE;
      break;

      case 'n':
      display_mode = DISPLAY;
      break;

      case 'g':
      grav *= 0.9;
      cout << "decrease gravity constant to "<<grav<<endl;
      cfd_demo->updateGrav(grav);
      break;

      case 'G':
      grav *= 1.0/0.9;
      cout << "increase gravity constant to "<<grav<<endl;
      cfd_demo->updateGrav(grav);
      break;

      case 'q':
      h *= 0.9;
      cout << "decreaese time step to "<<h<<endl;
      cfd_demo->updateTimeStep(h);
      break;

      case 'Q':
      h *= 1.0/0.9;
      cout << "increase time step to "<<h<<endl;
      cfd_demo->updateTimeStep(h);
      break;

      case 'b':
      has_boundary = !has_boundary;
      break;

      //case 'e':
      //writeTestImage( "test.jpg", density_map );
      //exit(0);
      //break;

      case ' ':
      if(is_on){
         is_on = false;
         cout << "stop simulation" << endl;
      }else{
         is_on = true;
         cout << "start simulation" << endl;
      }
      break;

      default:
      break;
   }
}

void cbMouseDown( int button, int state, int x, int y )
{
   if( button != GLUT_LEFT_BUTTON ) { return; }
   if( state != GLUT_DOWN ) { return; }
   xmouse_prev = x;
   ymouse_prev = y;
   DabSomePaint( x, y );
}

void cbMouseMove( int x, int y )
{
   xmouse_prev = x;
   ymouse_prev = y;
   DabSomePaint( x, y ); 
}

void PrintUsage()
{
   cout << "cfd_paint keyboard choices\n";
   cout << "space   start/stop simulation\n";
   cout << "b       turns on/off has_boundary\n";
   cout << "s       turns on painting source strength\n";
   cout << "o       turns on painting obstructions\n";
   cout << "+/-     increase/decrease brightness of display\n";
   cout << "r       resets brightness to default\n";
   cout << "c       toggles screen capture on/off\n";
   cout << "q       decrease size of the time step\n";
   cout << "Q       increase size of the time step\n";
   cout << "g       decrease gravity constant\n";
   cout << "G       increase gravity constant\n";
   cout << "e       quits the program\n";
}

//---------------------------------------------------

int main(int argc, char** argv)
{
    frame = 1;
    CmdLineFind clf( argc, argv );

    iwidth = clf.find("-iwidth", 512, "Horizontal width" );
    iheight = clf.find("-iheight", 512, "Vertical height" );

    nx = clf.find("-NX", 512, "Horizontal grid points" );
    ny = clf.find("-NY", nx, "Vertical grid points" );

    gs = clf.find("-GS", 30, "GS loops");
    iop = clf.find("-IOP", 5, "IOP loops");
        
    capture_screen = clf.findFlag("-capture");
    captured_file_basename = clf.find("-fname", "cfdbeginning" );

    setNbCores(4);
    string imagename = clf.find("-image", "grumpy.jpg", "Image to drive color");

    clf.usage("-h");
    clf.printFinds();
    PrintUsage();

    lx = float(iwidth);
    ly = float(iheight);
     
    is_on = false;
    size = nx * ny;
    scaling_factor = 1;

    paint_mode = PAINT_SOURCE;
    display_mode = DISPLAY;
    has_boundary = true;

    cfd_demo = new cfd(lx, ly, nx, ny, gs, iop, h, grav);

    display_map = new float[3*size];
    Initialize(display_map, 3*size, 0.0 );

    source_map = new float[size];
    Initialize(source_map, size, 0.0);

    obstruction_map = new float[size];
    Initialize(obstruction_map, size, 1.0);

    divergence_source_map = new float[size];
    Initialize(divergence_source_map, size, 0.0);

    if( imagename != "" )
    {
      readOIIOImage( imagename.c_str(), cfd_demo->getColorMap() );
    }

 
   InitializeBrushes();

   
   // GLUT routines
   glutInit(&argc, argv);

   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
   glutInitWindowSize( iwidth, iheight );

   // Open a window 
   char title[] = "cfd Demo";
   glutCreateWindow( title );
   
   glClearColor( 1,1,1,1 );

   glutDisplayFunc(&cbDisplay);
   glutIdleFunc(&cbIdle);
   
   glutKeyboardFunc(&cbOnKeyboard);
   glutMouseFunc( &cbMouseDown );
   glutMotionFunc( &cbMouseMove );
  
   glutMainLoop();
   return 1;
};
