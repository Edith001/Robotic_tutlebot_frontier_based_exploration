#include<cstdio>
#include<cstring>
#include<cctype>
#include<cstdlib>
#include<cmath>
#include<math.h>
#include<vector>
#include<iostream>
#include<string>
#include<queue>
#include<map>
#include<set>
#include<stack>
#include<algorithm>
#include<utility>
#include<sstream>
#include <random>
#include <fstream>

using namespace std ;

double nCmdsPerSecond =  10 ;

// ----------------- Random Generators --------------------

default_random_engine lgenerator(234254);
default_random_engine agenerator(412356);
normal_distribution<double> linearErrorModel(0.0,0.006) ;
normal_distribution<double> angularErrorModel(0.0,0.0002) ;

// ----------------- 2d  Vector Code Starts ----------------

#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)
#define EPS 1e-10
// +0.0000000001
#define S(x)	((x)*(x))
#define Z(x)	(fabs(x) < EPS)

struct vector2d{
	double x,y ;
	vector2d(){x = 0 ; y = 0; };
	vector2d(double _x , double _y ){ x = _x ; y=_y ;}
	bool isZero(){
		return Z(x) && Z(y);
	}
	double mag(){
		return sqrt(mag2());
	}
	double mag2(){
		return S(x)+S(y);
	}
};

vector2d operator-(vector2d a){ return vector2d(-a.x , -a.y );}
vector2d operator-(vector2d a ,  vector2d b){ return vector2d(a.x-b.x , a.y-b.y) ;}
vector2d operator+(vector2d a ,  vector2d b){ return vector2d(a.x+b.x , a.y+b.y) ;}
vector2d operator/(vector2d a , double b){ return vector2d(a.x/b , a.y/b);}
vector2d operator*(double a , vector2d b){ return vector2d(a*b.x,a*b.y);}
double operator*(vector2d a , vector2d b){ return a.x*b.y - a.y*b.x ;}
double operator^(vector2d a , vector2d b){ return a.x*b.x + a.y*b.y ;}

double mag( vector2d a) { return sqrt( S(a.x) + S(a.y)) ;}
double dist( vector2d a , vector2d b) {return sqrt( S(a.x-b.x) + S(a.y-b.y)) ;}

vector2d unit( vector2d a ){ return a/mag(a) ;}
double proj( vector2d a , vector2d b ){ return a^unit(b) ;}
vector2d projv( vector2d a , vector2d b ){ vector2d ub = unit(b) ; return (a^ub)*ub ;}
vector2d rotate( vector2d a , double ang ){ vector2d b = vector2d(-a.y,a.x) ; return cos(ang)*a+sin(ang)*b ;}


double dotp(vector2d a, vector2d b){
	return a.x*b.x + a.y*b.y ;
};

double crossp(vector2d a, vector2d b){
	return a.x*b.y - a.y*b.x;
};

double angle(vector2d a, vector2d b){
	//if(fabs(a.mag2())<EPS||fabs(b.mag2())<EPS)return 0;
	double v= dotp(a,b) / (a.mag()*b.mag()) ;
	v=MIN(v,1);
	v=MAX(v,-1);
	return acos(v);
};

//  Gives an angle in the range -pi to pi ( which basically covers everything)
double getAngle(vector2d a , vector2d b ){
	double v = crossp(a, b);
	if(v<0) return -1*angle(a,b);
	else return angle(a,b);
};

// ----------------- 2d  Vector Code Ends  ----------------


// ----------------- Particle Structure    ----------------

struct particle{
	vector2d location  ;
	vector2d direction ;
	double weight ;
};

particle moveParticle( particle par , double linearVel , double angularVel ){
	particle ret   ;
	vector2d vectorAdd = unit(par.direction) ;
	double rotAngle  = angularVel*(1 / nCmdsPerSecond ) + angularErrorModel(agenerator) ;
	double linearDisplacement = linearVel*(1 / nCmdsPerSecond ) + linearErrorModel(lgenerator) ;
	vectorAdd = linearDisplacement * rotate(vectorAdd,rotAngle) ;
	ret.location = par.location + vectorAdd ;
	ret.direction= rotate(par.direction ,rotAngle  ) ;
	return ret ;
}





// ----------------- Particle Structure    ----------------
#define particleN 20
particle presentParticle[particleN] ;
particle bufferParticle[particleN] ;

void printParticles( ofstream &outfile ){

	for( int i = 0 ; i< particleN ; i++ ){
		outfile<<presentParticle[i].location.x<<","<<presentParticle[i].location.y<<","<<presentParticle[i].direction.x<<","<<presentParticle[i].direction.y<<","<<endl ;
	}


}

int main(){


	ofstream outfile;
	ifstream infile ;
	double linearVel , angularVel ;
	//infile.open("squaretest.txt") ;
	infile.open("translate5m.txt") ;
	//outfile.open("particles.csv", std::ios_base::app);
	outfile.open("data.txt");
	int linecount = 0 ;

	for ( int i =0 ; i<particleN ; i++ ){
		presentParticle[i].direction = vector2d(1,0);
	}

	if( infile.is_open()){
		while( infile>>linearVel>>angularVel ){
			//cout<<linearVel<<angularVel<<endl ;
			for ( int i =0 ; i<particleN ; i++ ){
				presentParticle[i] = moveParticle(presentParticle[i] , linearVel , angularVel );

			}
			linecount++ ;

			if( linecount%20 == 0 ){
				printParticles(outfile);
			}
		}
	}

	outfile.close();
	infile.close();


	return 0 ;
}
