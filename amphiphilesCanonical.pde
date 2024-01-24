// AmphiphilesWorkshop
// by R M L Evans last updated 18/01/2024 (and initial coding 20/07/2010).
//
// Implements equilibrium Metropolis Monte Carlo algorithm for a model of amphiphiles on a grid.


float AA = -1.0;  // Hydrophobic-hydrophobic interation energy
float BB = 0.0;  // Hydrophilic-hydrophilic interaction energy
float AB = 0.0;   // Hydrophobic-hydrophilic interaction energy
                  // N.B. All interactions with the "solvent" (empty sites) have zero energy.

float phi = 0.1;   // Initial concentration
float kT = 0.2;    // Initial temperature

int Lx = 100;  // Grid size
int Ly = 60;

int[][] occ = new int[Lx][Ly];    // Occupancy array
int[][] direct = new int[Lx][Ly]; // Direction (orientation) array

int N=0;

//=============================================================================================
// Housekeeping for input fields:
import controlP5.*;
ControlP5 cp5;
String textValue = "";
void controlEvent(ControlEvent theEvent) {
  if(theEvent.isAssignableFrom(Textfield.class)){
    if (theEvent.getName() == "Concentration"){
      phi = float(theEvent.getStringValue());   
      updateNumberOfAmphiphiles();
    }
    if (theEvent.getName() == "Temperature"){
      kT = float(theEvent.getStringValue());
    }
  }
}

//=============================================================================================
// This "method" (i.e. function) is called once at start-up:

void setup() 
{
//  size(10*Lx,10*Ly+70);  // Variable window sizes are no longer allowed by Processing
  size(1000,670);  // Create the window.

// Set up the text input boxes:
  PFont font = createFont("arial",20); textFont(font);
  cp5 = new ControlP5(this); 
  cp5.addTextfield("Concentration").setPosition(5,10*Ly+18).setAutoClear(false);
  cp5.get(Textfield.class,"Concentration").setText(nf(phi,0,0)); 
  cp5.addTextfield("Temperature").setPosition(340,10*Ly+18).setAutoClear(false);     
  cp5.get(Textfield.class,"Temperature").setText(nf(kT,0,0)); 
  
  updateNumberOfAmphiphiles();
}

//---------------------------------------------------------------------------------------------
// This method is called repeatedly, by a timer:

void draw()
{  
  MonteCarloStep();
  updateGraphics();
}

//---------------------------------------------------------------------------------------------

void updateNumberOfAmphiphiles()
{
  int i, j, NewN;
  NewN = (int)(Lx*Ly*phi);
  while (NewN > N){
    i = (int) random(Lx);
    j = (int) random(Ly);
    if (occ[i][j]==0){  // If the randomly-chosen site is empty, occupy it with a randomly oriented particle
      occ[i][j] = 1;
      direct[i][j] = int(random(4));
      N++;
    }
  }
  while (NewN < N){
    i = (int) random(Lx);
    j = (int) random(Ly);
    if (occ[i][j]==1){  // If the randomly-chosen site is occupied, empty it
      occ[i][j] = 0;
      N--;
    }
  }
}
//---------------------------------------------------------------------------------------------

void updateGraphics(){
  background(160,180,255);
  
  for (int i=0; i<Lx; i++) for (int j=0; j<Ly; j++) if (occ[i][j] == 1)
  {
    fill(0);
    if (direct[i][j] == 0) triangle(i*10+10,j*10,i*10+10,j*10+10,i*10,j*10+10);
    else if (direct[i][j] == 1) triangle(i*10+10,j*10+10,i*10,j*10+10,i*10,j*10);
    else if (direct[i][j] == 2) triangle(i*10,j*10+10,i*10,j*10,i*10+10,j*10);
    else if (direct[i][j] == 3) triangle(i*10,j*10,i*10+10,j*10,i*10+10,j*10+10);
    fill(255);
    if (direct[i][j] == 0) triangle(i*10,j*10+10,i*10,j*10,i*10+10,j*10);
    else if (direct[i][j] == 1) triangle(i*10,j*10,i*10+10,j*10,i*10+10,j*10+10);
    else if (direct[i][j] == 2) triangle(i*10+10,j*10,i*10+10,j*10+10,i*10,j*10+10);
    else if (direct[i][j] == 3) triangle(i*10+10,j*10+10,i*10,j*10+10,i*10,j*10);
  }
}  

//---------------------------------------------------------------------------------------------

void MonteCarloStep()
{
  int i,j,i1,j1,i2,j2,d,oc,d1,oc1,d2,oc2;
  float dE,sum;

  for (int iteration=0; iteration<Lx*Ly; iteration++)
  {
// Pick a site:
    i = (int) random(Lx);
    j = (int) random(Ly);
// If it's occupied, try re-orienting it:
    if (occ[i][j]==1)
    {
      d = (direct[i][j] + int(random(3))) % 4; //Choose a new direction
      dE = Energy(i,j,1,d) - Energy(i,j,1,direct[i][j]);
      if ((dE < 0) || (random(1) < exp(-dE/kT)) ) direct[i][j] = d; // Accept new direction depending on Metropolis algorithm
    }

// Pick two non-contacting sites:
    i1 = (int) random(Lx);
    j1 = (int) random(Ly);
    do{
      i2 = (int) random(Lx);
      j2 = (int) random(Ly);
    } while (((abs(i2-i1)<2)||(abs(i2-i1)>(Lx-2)))&&((abs(j2-j1)<2)||(abs(j2-j1)>(Ly-2))));
// Try to swap them:
    oc1 = occ[i1][j1];
    oc2 = occ[i2][j2];
    d1 = direct[i1][j1];
    d2 = direct[i2][j2];
    dE = Energy(i1,j1,oc2,d2) + Energy(i2,j2,oc1,d1) - Energy(i1,j1,oc1,d1) - Energy(i2,j2,oc2,d2);
      if ((dE < 0) || (random(1) < exp(-dE/kT)) ){
        occ[i1][j1] = oc2;
        occ[i2][j2] = oc1;
        direct[i1][j1] = d2;
        direct[i2][j2] = d1;
      }
  }
} 

//---------------------------------------------------------------------------------------------
// Function (the model's Hamiltonian) to return the energy of all 4 nearest-neighbour
// interactions with site (i,j), if it is given occupancy oc (=0 or 1) and direction d (=0,1,2,3):
float Energy(int i,int j,int oc,int d)
{
  if (oc==1){
    // For each of the 4 neighbours, in the 4 directions relative to the orientation of the particle at site (i,j):
    int[] nA = new int[4];    // nA=1 if neighbour presents A end, nA=0 if it presents B or is empty
    int[] nB = new int[4];    // nB=1 if neighbour presents B end, nB=0 if it presents A or is empty
  
    nA[(4-d)%4] = occ[(i+1)%Lx][j] * (((1 + direct[(i+1)%Lx][j])%4)/2);
    nA[(5-d)%4] = occ[i][(j+1)%Ly] * (((4 + direct[i][(j+1)%Ly])%4)/2);
    nA[(6-d)%4] = occ[(i+Lx-1)%Lx][j] * (((3 + direct[(i+Lx-1)%Lx][j])%4)/2);
    nA[(7-d)%4] = occ[i][(j+Ly-1)%Ly] * (((2 + direct[i][(j+Ly-1)%Ly])%4)/2);
  
    nB[(4-d)%4] = occ[(i+1)%Lx][j] * (((3 + direct[(i+1)%Lx][j])%4)/2);
    nB[(5-d)%4] = occ[i][(j+1)%Ly] * (((2 + direct[i][(j+1)%Ly])%4)/2);
    nB[(6-d)%4] = occ[(i+Lx-1)%Lx][j] * (((1 + direct[(i+Lx-1)%Lx][j])%4)/2);
    nB[(7-d)%4] = occ[i][(j+Ly-1)%Ly] * (((4 + direct[i][(j+Ly-1)%Ly])%4)/2);
    
    return AA*(nA[0]+nA[1]) + AB*(nA[2]+nA[3]+nB[0]+nB[1]) + BB*(nB[2]+nB[3]);
  }
  else return 0.0;
}
