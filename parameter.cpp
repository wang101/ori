class para
{      public:
            double lambda;//meter 
            double pix_size;//perimeter of each pixel
            int demo_size; //space of 3D model
            int scr_size;//screen size
            double wave_number;//wavenumber = 1/lambda
            double distance;//distance from crystal to screen
            double epsilon;//parameter reserved for calculation of crystal
            char eulerdef[3]; //def of euler angle 'zxz' as default;
            int natom;//atom number
            bool magnetic_field;
           double magnetic_angle[3];
           double tesla;
           double miu;
           double temperature;
           para()
           {
  	eulerdef[0]='z';
            eulerdef[1]='x';
            eulerdef[2]='z';
            demo_size=15;
            scr_size=128;
            pix_size=0.008645052112070126;
            distance=1.0;
            lambda = 1.00e-10;
            natom = 3344;
             magnetic_field = false;
             magnetic_angle[0] = 0;
             magnetic_angle[1] = 0;
             magnetic_angle[2] = 1;
             tesla = 1.0;
             miu = -9.274e-21;
            temperature = 273;
            }};