#include "../../common/include/B4cDetectorConstruction.hh"
#include "../../common/include/B4cCalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4EllipticalTube.hh"



#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "math.h"

//Absorber Box    12.5''
//Path  /Users/jijunchen/Documents/Documents/Geant4_work/12P5_SSC/Step2/src


//读懂代码建议:不要先看变量，先看physical volume是什么(要结合可视图),physical volume是由什么logical volume定义的，
//logical volume是什么G4的变量定义的,然后具体看变量是怎么被定义的，一步一步往回看，之所以倒着看，是因为一开始的变量都是
//不可视的，变量的乱七八糟定义都是为了最后一步可视的的物理体，所以要先看见最终物理体是长什么样的。

//代码分成几个部分,如果只想看一部分，可以按照部分的分割线来注释掉，
//比如只看Top shielding,除了Top shielding和detector以外的部分都注释掉就可以了,
//注意！两个detector的部分是无论什么情况都不可以被注释掉,因为Eventaction里面的第71行,
//需要Get CollectionID("LaBrHitsCollection")和“HPGeHitsCollection



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4cDetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::B4cDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),
   fNofLayers(-1),
   fLaBrPV(nullptr),fHPGePV(nullptr)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::~B4cDetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::Construct()
{

  G4cout<<"Come here!!!!! Begining of the Construct~~~~~~~~~~~~~"<<G4endl;
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 

  G4Element *elH, *elD, *elHe, *elLi, *elBe, *elC, *elN;
  G4Element *elNe, *elAl, *elMn, *elFe, *elNi, *elCu, *elW, *elPb, *elU, *elO;
  G4Element *elK, *elHg, *elCa, *elNa, *elSi;

  G4String name,symbol;
  a=1.01*g/mole;
  elH=new G4Element(name="Hydrogen",symbol="H2",z=1.,a);

  a=2.01*g/mole;
  elD=new G4Element(name="Deuterium",symbol="D",z=1.,a);

  a=4.*g/mole;
  elHe=new G4Element(name="Helium",symbol="He",z=2.,a);

  a=6.94*g/mole;
  elLi=new G4Element(name="Lithium",symbol="Li",z=3.,a);

  a=9.01*g/mole;
  elBe=new G4Element(name="Beryllium",symbol="Be",z=4.,a);

  a=12.01*g/mole;
  elC=new G4Element(name="Carbon",symbol="C",z=6.,a);

  a=14.01*g/mole;
  elN=new G4Element(name="Nitrogen",symbol="N2",z=7.,a);

  a=16.*g/mole;
  elO=new G4Element(name="Oxygen",symbol="O2",z=8.,a);

  a=20.18*g/mole;
  elNe=new G4Element(name="Neon",symbol="Ne",z=10.,a);

  a=22.99*g/mole;
  elNa=new G4Element(name="Sodium",symbol="Na",z=11.,a);

  a=26.98*g/mole;
  elAl=new G4Element(name="Aluminium",symbol="Al",z=13.,a);

  a=28.085*g/mole;
  elSi=new G4Element(name="Silicon",symbol="Si",z=14.,a);

  a = 39.1*g/mole;
  elK  = new G4Element(name="K"  ,symbol="K" , z=19 , a);

  a=40.08*g/mole;
  elCa=new G4Element(name="Calcium",symbol="Ca",z=20.,a);

    
  a=54.94*g/mole;
  elMn=new G4Element(name="Manganese",symbol="Mn",z=25,a);
    
  a=55.850*g/mole;
  elFe=new G4Element(name="Iron",symbol="Fe",z=26.,a);
    
  a=58.7*g/mole;
  elNi=new G4Element(name="Nickel", symbol="Ni", z=28.,a);
    
  a=63.54*g/mole;
  elCu=new G4Element(name="Copper",symbol="Cu",z=29.,a);

  a=183.85*g/mole;
  elW=new G4Element(name="Tungsten",symbol="W",z=74.,a);

  a = 200.59*g/mole;
  elHg  = new G4Element(name="Hg"  ,symbol="Hg" , z=80, a);

  a=207.19*g/mole;
  elPb=new G4Element(name="Lead",symbol="Pb",z=82.,a);

  a=238.03*g/mole;
  elU=new G4Element(name="Uranium",symbol="U",z=92.,a);


  

    

  density = 2.5*g/cm3;
  G4Material*  ShieldingConcrete = new G4Material(name="ShieldingConcrete",density,10);
  

  new G4Material("Germanium", z=32., a= 72.64*g/mole, density= 5.323*g/cm3);
  new G4Material("Beryllium", z=4.,  a= 9.01*g/mole,  density= 1.850*g/cm3);
  new G4Material("Lead",      z=82., a= 207.2 * g/mole, density= 11.34*g/cm3);
  new G4Material("Copper",    z=29., a= 63.546* g/mole, density= 8.96 *g/cm3);
  new G4Material("Aluminum",  z=13., a= 26.98 * g/mole, density= 2.70 *g/cm3);
  new G4Material("Tungsten",  z=74., a= 183.84* g/mole, density= 19.25*g/cm3);
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
  
  ShieldingConcrete->AddElement(elH , 0.01);
  ShieldingConcrete->AddElement(elO , 0.529);
  ShieldingConcrete->AddElement(elNa , 0.016);
  ShieldingConcrete->AddElement(elHg , 0.002);
  ShieldingConcrete->AddElement(elAl , 0.034);
  ShieldingConcrete->AddElement(elSi , 0.337);
  ShieldingConcrete->AddElement(elK , 0.013);
  ShieldingConcrete->AddElement(elCa , 0.044);
  ShieldingConcrete->AddElement(elFe , 0.014);
  ShieldingConcrete->AddElement(elC , 0.001);

  /* //              97% Tugnsten
  density=18.5*g/cm3;
  G4Material*  Tungsten_97_per= new G4Material(name="Tungsten_97_per", density, 3);
  Tungsten_97_per->AddElement(elNi, 0.021);
  Tungsten_97_per->AddElement(elFe, 0.009);
  Tungsten_97_per->AddElement(elW,   0.97);
  */

  density=7.85*g/cm3;
  G4Material* Steel=new G4Material(name="Steel",density, 5);
  Steel->AddElement(elFe, 0.98);
  Steel->AddElement(elMn, 0.01);
  Steel->AddElement(elC, 0.004);
  Steel->AddElement(elCu, 0.003);
  Steel->AddElement(elSi, 0.003);
  
    
    

   G4Element* elLa = new G4Element("Lanthanum","La" ,57,138.905 *g/mole);
   G4Element* elBr = new G4Element("Bromine"  ,"Br" ,35, 79.904 *g/mole);
   G4Material* LaBr3 = new G4Material("LaBr3",5.08*g/cm3,2);
   LaBr3 -> AddElement(elLa,1);
   LaBr3 -> AddElement(elBr,3);

   G4Element* elCe  = new G4Element("Cerium"  ,"Ce" , z= 58.,a =140.116*g/mole);

   G4Material* LaBr3Ce = new G4Material("LaBr3Ce",5.08*g/cm3,2);
   LaBr3Ce -> AddElement(elCe,1.9*perCent);
   LaBr3Ce -> AddMaterial(LaBr3,98.1*perCent);
 

   G4Element* elB  = new G4Element("Boron"  ,"B" , 5., 10.81*g/mole);
   G4double ratio = 5;

   G4Material* polyethylene = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");

   G4Material* BP = new G4Material("BP",1.04*g/cm3,2);
   BP -> AddElement(elB,ratio*perCent);
   BP -> AddMaterial(polyethylene,(100-ratio)*perCent);



  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* B4cDetectorConstruction::DefineVolumes()
{
  // Geometry parameters

  G4double worldSizeXY =   50. *cm;
  G4double worldSizeZ  =   120. *cm; 
 
 
  auto LaBrMaterial   = G4Material::GetMaterial("LaBr3Ce");
  auto HPGeMaterial   = G4Material::GetMaterial("Germanium");

  auto WindowMaterial = G4Material::GetMaterial("Aluminum");
  auto EndcupMaterial = G4Material::GetMaterial("Aluminum");

  auto WorldMaterial    = G4Material::GetMaterial("Galactic");
  auto TungstenMaterial = G4Material::GetMaterial("Tungsten");
//auto TungstenMaterial = G4Material::GetMaterial("Tungsten_97_per");
    
  auto LeadMaterial     = G4Material::GetMaterial("Lead");
  auto CopperMaterial   = G4Material::GetMaterial("Copper");
  auto AluminumMaterial = G4Material::GetMaterial("Aluminum");
  auto BPMaterial = G4Material::GetMaterial("BP");
  auto ConcreteMaterial = G4Material::GetMaterial("ShieldingConcrete");
  auto AbsorberMaterial = G4Material::GetMaterial("G4_POLYETHYLENE");
  auto SteelMaterial    =G4Material::GetMaterial("Steel");
  //auto AbsorberMaterial = G4Material::GetMaterial("Air");

//////////////////////////////////////////////////////////////////////////////////////// 

  //     
  // World
  //
  G4Box* worldS 
    = new G4Box("World",           // its name
                 worldSizeXY*2.5, worldSizeXY*2.5, worldSizeZ*2.5); // its size
                         
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 WorldMaterial,  // its material
                 "World");         // its name
                                   
  G4PVPlacement* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


    
////////////////////////////////////////////////////////////////////////////////////////

  G4double offset_Y = 0;
  G4double offset_Spot = 40.6 * mm;
  G4double offset_HPGe = 2.45 * mm;  // beam hit the highest efficiceny of HPGe, JJ and Gage define
//G4double offset_HPGe = -21.6 * mm; // beam hit the left edge of crystal, JJ and Gage define

//G4double offset_HPGe = 12.0 * mm;  // beam hit the highest efficiceny of HPGe, Haichuan define
//G4double offset_HPGe = -20.0 * mm; // beam hit the left edge of crystal, Haichuan define
    
  G4double offset_LaBr = 00.0 * mm;


////////////////////////////////////////////////////////////////////////////////////////

  //     
  // LaBr Detector
  //

  G4double LaBr_EndcupR  = 4.8   *cm; //Radius
  G4double LaBr_EndcupL  = 6.63  *cm; //Length
  G4double LaBr_CrystalR = 3.81  *cm;
  G4double LaBr_CrystalL = 7.62  *cm;
  G4double Z_LaBr =  318*mm;

  G4double WindowD  = 1.5   *mm;
  G4double EndcupD  = 1.5   *mm;
  G4double AirD     = 4     *mm;


  G4double forward=1*25.4/2*mm; //forward doesn't have physics meaning, just a later change to move everything forward
  G4double gap = 0*mm;
//LaBr Crystal

  G4Tubs* LaBr_Detector = new G4Tubs("LaBr_Crystal", 0, LaBr_CrystalR, LaBr_CrystalL/2, 360.*deg, 360.*deg);
  G4LogicalVolume* fLaBrLV = new G4LogicalVolume(LaBr_Detector, LaBrMaterial, "fLaBrLV");
 
  fLaBrPV =
  new G4PVPlacement(
                 0,
                 G4ThreeVector(offset_Spot+offset_LaBr, 0., Z_LaBr+WindowD+AirD+LaBr_CrystalL/2 - forward + gap),
                 //G4ThreeVector(2.07*25.4*mm,0., 13.6*25.4*mm), //from engineer data
                  fLaBrLV,
                 "fLaBrPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


//LaBr EndCup

  G4Tubs* LaBr_Endcup = new G4Tubs("LaBr_Endcup", LaBr_EndcupR-EndcupD, LaBr_EndcupR, LaBr_EndcupL/2, 360.*deg, 360.*deg);
  G4LogicalVolume* LaBr_fEndcupLV = new G4LogicalVolume(LaBr_Endcup, EndcupMaterial, "LaBr_fEndcupLV");

  G4PVPlacement* LaBr_fEndcupPV =
  new G4PVPlacement(
                 0,
                 G4ThreeVector(offset_Spot+offset_LaBr, 0., Z_LaBr+LaBr_EndcupL/2 - forward + gap),
                 LaBr_fEndcupLV,
                 "LaBr_fEndcupPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


//LaBr Window

  G4Tubs* WindowLaBr   = new G4Tubs("WindowLaBr", 0, LaBr_EndcupR-EndcupD, WindowD/2, 360.*deg, 360.*deg);
  G4LogicalVolume* WindowLVLaBr = new G4LogicalVolume(WindowLaBr, WindowMaterial, "WindowLVLaBr");

  G4PVPlacement* fWindowPVLaBr =
  new G4PVPlacement(
                 0,
                 G4ThreeVector(offset_Spot+offset_LaBr, 0., Z_LaBr+WindowD/2-forward + gap),
                 WindowLVLaBr,
                 "fWindowPVLaBr",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);



////////////////////////////////////////////////////////////////////////////////////////

  //     
  // HPGe Detector
  //

  G4double HPGe_EndcupR  = 3.95*25.4/2*mm;
  G4double HPGe_EndcupL  = 173  *mm;
  G4double HPGe_CrystalR = 72.1/2*mm;
  G4double HPGe_CrystalL = 78.5  *mm;
//G4double Z_HPGe =  98.61*mm;  // For Edge
  G4double Z_HPGe =  63.131*mm; // For High efficiency
    
    
    
    
    
  G4double HoleR    = 10.5/2*mm;
  G4double HoleL    = 64.7  *mm;

  G4double Capsule_Wallthick = 0.8*mm;
  G4double Capsule_Windowthick = 0.03*mm;
  G4double Capsule_Endthick = 3*mm;
  G4double Capsule_Walllength = 130*mm;


  G4RotationMatrix* rotBox1 = new G4RotationMatrix();
  G4RotationMatrix* rotBox2 = new G4RotationMatrix();
  rotBox1->rotateY(45*deg);
  rotBox2->rotateX(90*deg);


//Crystal
  G4Tubs* HPGe_Hole1 = new G4Tubs  ("HPGe_Hole1", 0, HoleR, HoleL/2, 360.*deg, 360.*deg);
  G4Sphere* HPGe_Hole2 = new G4Sphere("HPGe_Hole2", 0, HoleR, 0*deg, 180*deg, 0*deg, 180*deg);
  G4UnionSolid* HPGe_Hole = new G4UnionSolid("HPGe_Hole", HPGe_Hole1, HPGe_Hole2, rotBox2, G4ThreeVector(0, 0, -HoleL/2));

  G4Tubs* HPGe_Crystal = new G4Tubs("HPGe_Crystal", 0, HPGe_CrystalR, HPGe_CrystalL/2, 360.*deg, 360.*deg);
  G4SubtractionSolid* HPGe_Detector = new G4SubtractionSolid("HPGe_Detector", HPGe_Crystal, HPGe_Hole, 0, G4ThreeVector(0, 0, (HPGe_CrystalL-HoleL)/2));


  G4LogicalVolume* fHPGeLV = new G4LogicalVolume(HPGe_Detector, HPGeMaterial, "fHPGeLV");

  fHPGePV =
  new G4PVPlacement(
                 rotBox1,
                 G4ThreeVector(-offset_Spot + offset_HPGe - (WindowD + AirD + Capsule_Windowthick + HPGe_CrystalL/2)*sqrt(2)/2, 0., Z_HPGe + (WindowD + AirD + Capsule_Windowthick + HPGe_CrystalL/2)*sqrt(2)/2),
                 fHPGeLV,
                 "fHPGePV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


  G4LogicalVolume* HoleLV = new G4LogicalVolume(HPGe_Hole, WorldMaterial, "HoleLV");

  HolePV =
  new G4PVPlacement(
                 rotBox1,
                 G4ThreeVector(-offset_Spot + offset_HPGe - (WindowD + AirD + Capsule_Windowthick + HPGe_CrystalL - HoleL/2)*sqrt(2)/2, 0., Z_HPGe + (WindowD + AirD + Capsule_Windowthick + HPGe_CrystalL - HoleL/2)*sqrt(2)/2),
                 HoleLV,
                 "HolePV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

  //EndCup
  G4Tubs* HPGe_Endcup = new G4Tubs("HPGe_Endcup", HPGe_EndcupR-EndcupD, HPGe_EndcupR, HPGe_EndcupL/2, 360.*deg, 360.*deg);
  G4LogicalVolume* HPGe_fEndcupLV = new G4LogicalVolume(HPGe_Endcup, EndcupMaterial, "HPGe_fEndcupLV");

  G4PVPlacement* HPGe_fEndcupPV =
  new G4PVPlacement(
                 rotBox1,
                 G4ThreeVector(-offset_Spot + offset_HPGe - HPGe_EndcupL*sqrt(2)/4, 0., Z_HPGe + HPGe_EndcupL*sqrt(2)/4),
                 HPGe_fEndcupLV,
                 "HPGe_fEndcupPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


//Window
  G4Tubs* WindowHPGe   = new G4Tubs("WindowHPGe",0, HPGe_EndcupR-EndcupD, WindowD/2, 360.*deg, 360.*deg);
  G4LogicalVolume* WindowLVHPGe = new G4LogicalVolume(WindowHPGe,WindowMaterial,"WindowLVHPGe");


  G4PVPlacement* fWindowPVHPGe =
  new G4PVPlacement(
                 rotBox1,
                 G4ThreeVector(-offset_Spot + offset_HPGe, 0.,Z_HPGe+WindowD/2),
                 //G4ThreeVector(-offset_Spot + offset_HPGe, 0.,Z_HPGe+WindowD/2),
                 WindowLVHPGe,
                 "fWindowPVHPGe",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


//Capsule
  G4Tubs* HPGe_Capsule1 = new G4Tubs("HPGe_Capsule1", HPGe_CrystalR, HPGe_CrystalR+Capsule_Wallthick, Capsule_Walllength/2, 360.*deg, 360.*deg);
  G4Tubs* HPGe_Capsule2 = new G4Tubs("HPGe_Capsule2", 0, HPGe_CrystalR, Capsule_Windowthick/2, 360.*deg, 360.*deg);
  G4Tubs* HPGe_Capsule3 = new G4Tubs("HPGe_Capsule3", 0, HPGe_CrystalR, Capsule_Endthick/2, 360.*deg, 360.*deg);
  G4Tubs* HPGe_Capsule4 = new G4Tubs("HPGe_Capsule4", 0, HoleR, Capsule_Endthick/2 + 0.001, 360.*deg, 360.*deg);

  G4UnionSolid* HPGe_Capsule12 = new G4UnionSolid("HPGe_Capsule12", HPGe_Capsule1, HPGe_Capsule2, 0, G4ThreeVector(0, 0, -(Capsule_Walllength - Capsule_Windowthick)/2));
  G4SubtractionSolid* HPGe_Capsule34 = new G4SubtractionSolid("HPGe_Capsule34", HPGe_Capsule3, HPGe_Capsule4, 0, G4ThreeVector(0, 0, 0));

  G4LogicalVolume* HPGe_CapsuleLV12 = new G4LogicalVolume(HPGe_Capsule12, AluminumMaterial, "HPGe_CapsuleLV12");
  G4LogicalVolume* HPGe_CapsuleLV34 = new G4LogicalVolume(HPGe_Capsule34, AluminumMaterial, "HPGe_CapsuleLV34");

    G4PVPlacement* HPGe_CapsulePV12 =
    new G4PVPlacement(
                 rotBox1,
                 G4ThreeVector(-offset_Spot + offset_HPGe - (WindowD + AirD + Capsule_Walllength/2)*sqrt(2)/2, 0., Z_HPGe + (WindowD + AirD + Capsule_Walllength/2)*sqrt(2)/2),
                 HPGe_CapsuleLV12,
                 "HPGe_CapsulePV12",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


    G4PVPlacement* HPGe_CapsulePV34 =
    new G4PVPlacement(
                 rotBox1,
                 G4ThreeVector(-offset_Spot + offset_HPGe - (WindowD + AirD + Capsule_Walllength - Capsule_Endthick/2)*sqrt(2)/2, 0., Z_HPGe + (WindowD + AirD + Capsule_Walllength - Capsule_Endthick/2)*sqrt(2)/2),
                 HPGe_CapsuleLV34,
                 "HPGe_CapsulePV34",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


G4double height = 6.0*2.54*cm;

G4int Detectoroption=1; // cross-talk sheilding option


if(Detectoroption == 1) // cross-talk sheilding option
{

  G4double InnerFront_l = 4.7*25.4*mm; //4.7’‘ is 119.5mm   l mean length, but not necessary in x direction
  G4double InnerSide_l = 159.24*mm - forward; //forward doesn't have physics meaning, just a later change to move everything forward
  
 //G4double InnerSide_Z = 4.4*25.4*mm; //which is 111.76*mm;
   G4double InnerSide_Z = 5.42*25.4*mm;
    
  G4double Holesize= 1.57*25.4*mm/2;

  G4double  fHalf_Pb = 25*mm; //f meaning???? thickness of Pb?
  G4double  rHalf_Pb = (0.787*25.4)/2*mm;
  G4double  Half_Cu = (0.16*25.4)/2*mm;//
  G4double  Half_Cu_Mid_side= (0.16*25.4)/2*mm;
  G4double  Half_Al = (0.04*25.4)/2*mm;
  G4double  Half_Al_Mid_side = (0.04*25.4)/2*mm;
  G4double  thickHalf_Pb = 138.84*mm/2;


 //G4double distance = 159.24*mm + 2*fHalf_Pb - forward;
  G4double distance = 2*fHalf_Pb + 2.25*25.4*mm;
  G4double FrontAl_Z1  = InnerSide_Z + Half_Al;
  G4double FrontCu_Z1  = FrontAl_Z1 + Half_Al + Half_Cu;
  G4double FrontPb_Z1  = FrontCu_Z1 + Half_Cu + fHalf_Pb;

  G4double FrontPb_Z2  = FrontPb_Z1 + distance;
  G4double FrontCu_Z2  = FrontPb_Z2 + Half_Cu + fHalf_Pb;
  G4double FrontAl_Z2  = FrontCu_Z2 + Half_Al + Half_Cu;

  G4double SidePb_X  = -0.287*25.4*mm;//
  G4double SideCu_X1 = SidePb_X  + Half_Cu + rHalf_Pb;
  G4double SideCu_X2 = SidePb_X  - Half_Cu - rHalf_Pb;
  G4double SideAl_X1 = SideCu_X1 + Half_Cu + Half_Al;
  G4double SideAl_X2 = SideCu_X2 - Half_Cu - Half_Al;

  G4double Front_extX  =  16.384*cm;
  G4double InnerSideBL_Z= (16.94+3)*25.4*mm ; // Need number for gage
  G4double InnerSideBR_Z= 11.48*25.4*mm;
    
  G4double Side_extZ = 5.5*cm;

  G4double Finner_X = 6.4385*cm;
  G4double Finner_offset = 3.372*cm;

  G4double Side_moveZ = (fHalf_Pb + Half_Cu )*2;
  G4double InnerSideSExtZ = Z_LaBr - forward + InnerSideBL_Z/2;
  G4double InnerFrontExtX = Finner_X - (InnerFront_l+Front_extX)/2;

  G4double InnerSideSExtX_Pb = InnerFrontExtX - Front_extX/2 + thickHalf_Pb;
    
  G4double rTheta_Pb_1  =(-45)*degree;
  G4double rTheta_Pb_2  =(90)*degree;
  G4double thick_Pb   =5.26*25.4*mm; //z方向上的长度，
  G4double length_Pb_1=16.65*25.4*mm; // 16.94''-0.2828''(thickness of Cu and Al)
  G4double length_Pb_2=11.48*25.4*mm;
    
    
  G4double InnerMidCu_X       = SidePb_X  - rHalf_Pb - Half_Cu;
  G4double InnerBackLead_X    = SidePb_X  - rHalf_Pb - Half_Cu*2 - thick_Pb/2;
  G4double InnerSideBackRCu_X = SidePb_X  - rHalf_Pb - Half_Cu*2 - thick_Pb - Half_Cu;
  G4double InnerSideS_Z = InnerSide_Z + distance/2 +fHalf_Pb +2*Half_Cu ;
    
    
  G4double Relative_Hole     = -23.86*mm; //Need to check with gage
////////////////////////////////////////////////////////////////////////////////////////
  
  G4Box* LeadInnerS_F  = new G4Box("LeadInnerS_F", rHalf_Pb, height/2, (distance+fHalf_Pb*2)/2);
      

//G4Box* LeadInnerS_B  = new G4Box("LeadInnerS_B", thickHalf_Pb, height/2, InnerSideB_Z/2);
      
  G4Trap* LeadInnerS_B = new G4Trap( "LeadInnerS_B",
                                      height,  thick_Pb,
                                      length_Pb_1, length_Pb_2);
    
  G4RotationMatrix* rotLeadBlock = new G4RotationMatrix();
  rotLeadBlock->rotateX(90*deg);
  rotLeadBlock->rotateZ(-90*deg);

    
  G4Box* LeadInnerFLayer = new G4Box ("LeadInnerFLayer", InnerFront_l/2, height/2, fHalf_Pb);
  G4Tubs* LeadInnerFHole = new G4Tubs("LeadInnerFHole",  0, Holesize, fHalf_Pb*1.01, 360.*deg, 360.*deg); //All *1.01 in this code means to visualize the internal structure, in this case, it is the internal hole
  G4SubtractionSolid* LeadInnerF = new G4SubtractionSolid("LeadInnerF", LeadInnerFLayer, LeadInnerFHole, 0, G4ThreeVector( -23.86*mm, -offset_Y, 0));


  G4Box* LeadInnerF_ext  = new G4Box("LeadInnerF_ext", Front_extX/2, height/2, fHalf_Pb);


  G4LogicalVolume* LeadInnerSLV1 = new G4LogicalVolume(LeadInnerS_F, LeadMaterial, "LeadInnerSLV1");
  G4LogicalVolume* LeadInnerSLV2 = new G4LogicalVolume(LeadInnerS_B, LeadMaterial, "LeadInnerSLV2");



  G4LogicalVolume* LeadInnerFLV = new G4LogicalVolume(LeadInnerF, LeadMaterial, "LeadInnerFLV");
  G4LogicalVolume* LeadInnerF_extLV = new G4LogicalVolume(LeadInnerF_ext, LeadMaterial, "LeadInnerF_extLV");


      LeadInnerSPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(SidePb_X, offset_Y, InnerSideS_Z),
                 LeadInnerSLV1,
                 "LeadInnerSPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      LeadInnerSPV2 =
      new G4PVPlacement(
                 rotLeadBlock,
               //G4ThreeVector(0, 0, 0),
                 G4ThreeVector(SidePb_X-rHalf_Pb-Half_Cu*2-thick_Pb/2, offset_Y, 17.902*25.4*mm),
               //G4ThreeVector(-3.4*25.4*mm, offset_Y, 17.902*25.4*mm),//(InnerSideSExtX_Pb, offset_Y, InnerSideSExtZ - 100*mm)Need update!
                 LeadInnerSLV2,
                 "LeadInnerSPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


      G4PVPlacement* LeadInnerFPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(Finner_X, offset_Y, FrontPb_Z1),
                 LeadInnerFLV,
                 "LeadInnerFPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      G4PVPlacement* LeadInnerFPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(Finner_X, offset_Y, FrontPb_Z2),
                 LeadInnerFLV,
                 "LeadInnerFPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

////////////////////////////////////////////////////////////////////////////////////////

  G4Box* CopperInnerS_F  = new G4Box("CopperInnerS_F", Half_Cu, height/2, (distance+fHalf_Pb*2+Half_Cu*2)/2);
  G4Box* CopperInnerS_F_edge = new G4Box("CopperInnerS_F_edge", rHalf_Pb, height/2, Half_Cu);
    
    
  G4Box* CopperInnerS_BackLeft  = new G4Box("CopperInnerS_BackLeft", Half_Cu_Mid_side, height/2, InnerSideBL_Z/2);
  G4Box* CopperInnerS_BackRight = new G4Box("CopperInnerS_BackRight",Half_Cu_Mid_side, height/2, InnerSideBR_Z/2);

  G4Box* CopperInnerFLayer = new G4Box ("CopperInnerFLayer", InnerFront_l/2, height/2, Half_Cu);
  G4Tubs* CopperInnerFHole = new G4Tubs("CopperInnerFHole",  0, Holesize, Half_Cu*1.01, 360.*deg, 360.*deg);
  G4SubtractionSolid* CopperInnerF = new G4SubtractionSolid("CopperInnerF", CopperInnerFLayer, CopperInnerFHole, 0, G4ThreeVector( -23.86*mm, -offset_Y, 0));

    
    
  G4LogicalVolume* CopperInnerSLV0 = new G4LogicalVolume(CopperInnerS_F_edge, CopperMaterial, "CopperInnerSLV0");
  G4LogicalVolume* CopperInnerSLV1 = new G4LogicalVolume(CopperInnerS_F, CopperMaterial, "CopperInnerSLV1");
  G4LogicalVolume* CopperInnerSLV2 = new G4LogicalVolume(CopperInnerS_BackLeft, CopperMaterial, "CopperInnerSLV2");
  G4LogicalVolume* CopperInnerSLV3 = new G4LogicalVolume(CopperInnerS_BackRight, CopperMaterial, "CopperInnerSLV3");
    
  G4LogicalVolume* CopperInnerFLV = new G4LogicalVolume(CopperInnerF, CopperMaterial, "CopperInnerFLV");

    
    
    
    CopperInnerSPV0 =
    new G4PVPlacement(
               0,
               G4ThreeVector(SidePb_X, offset_Y, InnerSideS_Z-(distance+fHalf_Pb*2)/2-Half_Cu/2),
               CopperInnerSLV0,
               "CopperInnerSPV0",
               worldLV,
               false,
               0,
               fCheckOverlaps);
    
    
    
    

      CopperInnerSPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(SideCu_X1, offset_Y, InnerSideS_Z -2*mm ),
                 CopperInnerSLV1,
                 "CopperInnerSPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

    

      CopperInnerSPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(InnerMidCu_X, offset_Y, InnerSideSExtZ - 175*mm ),
                 CopperInnerSLV2,
                 "CopperInnerSPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      CopperInnerSPV3 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(InnerSideBackRCu_X, offset_Y, 19.1*25.4*mm ),
                 CopperInnerSLV3,
                 "CopperInnerSPV3",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      G4PVPlacement* CopperInnerFPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(Finner_X, offset_Y, FrontCu_Z1),
                 CopperInnerFLV,
                 "CopperInnerFPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      G4PVPlacement* CopperInnerFPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(Finner_X, offset_Y, FrontCu_Z2),
                 CopperInnerFLV,
                 "CopperInnerFPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);



//////////////////////////////////////////////////////////////////////////////////////// Sandwich shape ( rotation cut )

    
  G4double Block_Theta= atan(0.5);
  G4double Side_Theta = 45*degree;

  G4double Block_change=10*cm;
  G4double Block_length = InnerSide_l + Side_moveZ - Block_change; // Block_change means????
  //G4double BlockSide_Length =(InnerSideBL_Z-InnerSideBL_Z)*sqrt(2)*1.05; //The reason for *1.05 is to make sure the Sandwich shape can touch the other size of block seamlessly, it overlapped with other object, but since it is the copper overlapping the copper, it doesn't change any thing. But if it is the lead overlapping the copper, the overlapping part would be considered as the lead.
  G4double BlockSide_Length =7.32*25.4*mm;
    
 
  G4double Block_X = SidePb_X - rHalf_Pb - Half_Cu_Mid_side*2  - thick_Pb/2 - Half_Cu_Mid_side/sqrt(2);
    
  G4double Block_Z = 10.8*25.4*mm;

  G4RotationMatrix* rotBlock = new G4RotationMatrix();
  rotBlock->rotateY(90*deg);

  G4RotationMatrix* rotSide = new G4RotationMatrix();
  rotSide->rotateY(45*deg);


  G4Box* LeadBlockWhole = new G4Box("LeadBlockWhole", Block_length/2, Block_length/2, Block_length/2);
  G4Box* LeadBlockCut = new G4Box("LeadBlockCut", Block_length/sqrt(2), Block_length/sqrt(2), Block_length/sqrt(2));
  G4SubtractionSolid* LeadBlockShape = new G4SubtractionSolid("LeadBlockShape", LeadBlockWhole, LeadBlockCut, rotSide, G4ThreeVector(-Block_length/2, 0, Block_length/2));

  G4LogicalVolume* LeadBackBlockLV = new G4LogicalVolume(LeadBlockShape, LeadMaterial, "LeadBackBlockLV");

    
 G4Box* CopperBlockShape = new G4Box("CopperBlockShape", Half_Cu_Mid_side + 1, height/2, BlockSide_Length/2 + 1.5);
 G4LogicalVolume* CopperBackBlockLV = new G4LogicalVolume(CopperBlockShape, CopperMaterial, "CopperBackBlockLV");

 CopperBackBlockPV =
  new G4PVPlacement(
                 rotSide,
                 G4ThreeVector(Block_X , offset_Y, Block_Z),
                  CopperBackBlockLV,
                 "CopperBackBlockPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


////////////////////////////////////////////////////////////////////////////////////////

 G4double FrontEdge2_l = InnerSide_l + Side_moveZ - Block_length;
 G4Box* CopperFrontEdge2 = new G4Box("CopperFrontEdge2", Half_Cu, height/2, FrontEdge2_l/2);
 G4LogicalVolume* CopperFrontEdgeLV2 = new G4LogicalVolume(CopperFrontEdge2, CopperMaterial, "CopperFrontEdgeLV2");

 CopperFrontEdgePV2 =
  new G4PVPlacement(
                 0,
                 G4ThreeVector(SideCu_X2, offset_Y, InnerSide_Z + FrontEdge2_l/2),
                 CopperFrontEdgeLV2,
                 "CopperFrontEdgePV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

/*
 G4Box* AluminumFrontEdge2 = new G4Box("AluminumFrontEdge2", Half_Al, height/2, FrontEdge2_l/2);
 G4LogicalVolume* AluminumFrontEdgeLV2 = new G4LogicalVolume(AluminumFrontEdge2, AluminumMaterial, "AluminumFrontEdgeLV2");

 AluminumFrontEdgePV2 =
  new G4PVPlacement(
                 0,
                 G4ThreeVector(SideAl_X2, offset_Y, InnerSide_Z + FrontEdge2_l/2),
                 AluminumFrontEdgeLV2,
                 "AluminumFrontEdgePV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

////////////////////////////////////////////////////////////////////////////////////////
  G4double FrontEdge3_l = Front_extX -2*rHalf_Pb - FrontEdge2_l;


 G4Box* CopperFrontEdge3 = new G4Box("CopperFrontEdge3", FrontEdge3_l/2, height/2, Half_Cu);
 G4LogicalVolume* CopperFrontEdgeLV3 = new G4LogicalVolume(CopperFrontEdge3, CopperMaterial, "CopperFrontEdgeLV3");
 
  CopperFrontEdgePV3 =
  new G4PVPlacement(
                 0,
                 G4ThreeVector(SidePb_X - rHalf_Pb - FrontEdge2_l - FrontEdge3_l/2, offset_Y, InnerSide_Z + (InnerSide_l + Side_moveZ) - Half_Cu),
                 CopperFrontEdgeLV3,
                 "CopperFrontEdgePV3",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

 G4Box* AluminumFrontEdge3 = new G4Box("AluminumFrontEdge3", FrontEdge3_l/2, height/2, Half_Al);
 G4LogicalVolume* AluminumFrontEdgeLV3 = new G4LogicalVolume(AluminumFrontEdge3, AluminumMaterial, "AluminumFrontEdgeLV3");
 AluminumFrontEdgePV3 =
  new G4PVPlacement(
                 0,
                 G4ThreeVector(SidePb_X - rHalf_Pb - FrontEdge2_l - FrontEdge3_l/2, offset_Y, InnerSide_Z + (InnerSide_l + Side_moveZ) - Half_Cu*2 -Half_Al),
                 AluminumFrontEdgeLV3,
                 "AluminumFrontEdgePV3",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);
*/
}



//////////////////////////////////////////////////////

 
    

  G4double absorb_hW = 50/2*mm;   // Half of the width! Half
  G4double absorb_hH = 50/2*mm;
  G4double absorb_hT = 390/2*mm;
  G4double AbsorberZ = -2.5*m;
  G4double Whole_front_thickness_in_Z= 12 *25.4*mm; // -2 for 4 '' thinkness Tungsten and +2 for 8'' thickness Tungsten
  G4double Gap_Absorber_SSC= 0.2*25.4*mm;
 
    
  G4Box* Absorber = new G4Box("Absorber", absorb_hW, absorb_hH, absorb_hT);
   

  //G4double absorb_hW = 215.9*mm;
  //G4double absorb_hH = 228.6*mm;
  //G4double absorb_hT = 390/2*mm;
  //G4double AbsorberZ = -2.5*m;
  //G4double Absorber_Hole_r = 20*mm;
  //G4ThreeVector Absorber_offset = G4ThreeVector(offset_Spot, 0, 0);
  //G4EllipticalTube* Absorber_Hole = new G4EllipticalTube("Absorber_Hole", Absorber_Hole_r, Absorber_Hole_r, absorb_hT*1.02);
  //G4SubtractionSolid* Absorber = new G4SubtractionSolid("Absorber", Absorber_Shape, Absorber_Hole, 0, Absorber_offset);
    
  G4LogicalVolume* AbsorberLV = new G4LogicalVolume(Absorber, AbsorberMaterial, "AbsorberLV");


  AbsorberPV =
    new G4PVPlacement(
                 0,
                 G4ThreeVector(-offset_Spot, 0.,  - Whole_front_thickness_in_Z- Gap_Absorber_SSC - absorb_hT),
                 AbsorberLV,
                 "AbsorberPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

    
/*
//////////////////////////////////////////////////////

  G4double VD101_hW = 215.9*mm;
  G4double VD101_hH = 228.6*mm;
  G4double VD101_hT=1*mm;

    G4EllipticalTube* VD101 = new G4EllipticalTube("VD101", VD101_hW, VD101_hH, VD101_hT);
    G4LogicalVolume* VD101LV = new G4LogicalVolume(VD101, WorldMaterial, "VD101LV");

    VD101PV =
    new G4PVPlacement(
                 0,
                 G4ThreeVector(0., 0., -320.*mm),
                 VD101LV,
                 "VD101PV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);
*/

//////////////////////////////////////////////////////

G4double leaddepth     = 2*2.54*cm;
G4double leaddepth2    = 4*2.54*cm;
G4double copperdepth   = 0.5 * 2.54*cm;
G4double aluminumdepth = 0.1 * 2.54*cm;
G4double BPdepth       = 0.5 * 2.54*cm;

G4double thicker       = 2.5 * 2.54*cm;
G4double Wdepth1       = 3.6 * 2.54*cm + thicker;
G4double Bottom_thick  = aluminumdepth + copperdepth + leaddepth*2 + BPdepth*2;
 
   
//#########################Front Shielding Wall ###########################

  G4double fTheta=-atan(0.5);
  G4double front_length = 18.7*25.4*mm;
    
//#########################Tungsten Paramter#########################
  G4double delta_WLR_Plus = 0*25.4*mm; // -1.5*25.4*mm for 11 collimator, 0 for 12.5collimator , 1.5*25.4*mm for 14 collimator
  G4double delta_WlR      = 2.5*25.4*mm + delta_WLR_Plus; //All the delta_ means Sudden/Temporary change
  G4double delta_WlL      = 1.3*25.4*mm; //All the delta_ means Sudden/Temporary change
  G4double fW_x           =  -0.15*25.4*mm - delta_WlR/2 + delta_WlL/2;
//G4double fW_x           =  -0.75*25.4*mm;//position of tungsten in x direction for 12.5 tungsten is 0.75''
  G4double W_length       = 8.7*25.4*mm + delta_WlR + delta_WlL;
  G4double W_height       = 7*25.4*mm;
  G4double leak           = 0.4 * 25.4*mm;
//G4double Wdepth_f       = 2*25.4*mm; // f means front part of Tungsten
  G4double Wdepth_f       = 2*25.4*mm;
  G4double Wdepth_b       = 2*25.4*mm; // b means back part of Tungsten
  G4double Wdepth         = Wdepth_f + Wdepth_b ;
  G4double fPb_length     = 4.75*25.4*mm;
  G4double fPb_lengthL    = 4.75*25.4*mm; //Left side thickness in x direction is 4inche
  G4double fPb_lengthR    = fPb_lengthL;
  G4double fPb_lengthT    = 22.78*25.4*mm + delta_WLR_Plus;

   
    

  G4double fAl_x  = fPb_lengthT;
//G4double fAl_x  = 22.78*25.4*mm;  //unentangle with fPb_length
  G4double fCu_x  = fAl_x;
  G4double fBP1_x = fCu_x;
  G4double fPb_x  = fBP1_x;
  G4double fBP2_x = fPb_x;
  
  G4double f_delta_y = 4.93*25.4*mm;
    
  G4double front_dY  = 0.27*25.4*mm;
  G4double fAl_y     = 17.13*25.4*mm;
  G4double fCu_y     = fAl_y;
  G4double fBP1_y    = fCu_y;
  G4double fPb_y     = fBP1_y;
  G4double fBP2_y    = fPb_y;
    
  G4double Aperture_HPGe1 = 1.5*cm*cm;  // The Aperture Size for the HPGe Detector 前面的孔面积
  G4double Aperture_LaBr1 = 1.5*cm*cm;  // The Aperture Size for the LaBr Detector 前面的孔面积
  G4double Aperture_HPGe2 = 0.5*cm*cm;  // The Aperture Size for the HPGe Detector 后面的孔面积
  G4double Aperture_LaBr2 = 0.5*cm*cm;  // The Aperture Size for the LaBr Detector 后面的孔面积
  G4double Pi  = 3.14159265358979323846;

  //因为有collimator有plug所以分别有 r_HPGe 和 r_HPGe2
    
  G4double r_HPGe = sqrt(Aperture_HPGe1/Pi);
  G4double r_LaBr = sqrt(Aperture_LaBr1/Pi);
 
  G4double r_HPGe2 = sqrt(Aperture_HPGe2/Pi);
  G4double r_LaBr2 = sqrt(Aperture_LaBr2/Pi);

  //G4double r_HPGe2 = 0.36/2*25.4*mm;   //collimator后面的空半径
  //G4double r_LaBr2 = 0.55/2*25.4*mm;
    
  G4ThreeVector Pos_hole1 = G4ThreeVector( -offset_Spot - fW_x, -front_dY , 0);
  G4ThreeVector Pos_hole2 = G4ThreeVector( +offset_Spot - fW_x, -front_dY, 0);
  //  G4ThreeVector Pos_hole1 = G4ThreeVector((front_length + BPdepth*2 + leaddepth + copperdepth)/2 + aluminumdepth/4 - 8.8*25.4 - 81.2, -front_dY, 0);
  //  G4ThreeVector Pos_hole2 = G4ThreeVector((front_length + BPdepth*2 + leaddepth + copperdepth)/2 + aluminumdepth/4 - 8.8*25.4, -front_dY, 0);
    
  //以上两行代表挖掉了两个洞，一个是给LaBr的一个是给HPGe的洞洞
  //Spot_FXXX_HPGe 指的是collimator后面的铅块和铜块和铝块的孔径的大小,为了保证后面的大洞包住前面collimator的大小洞，所以后面的大洞
  //的r_HPGe和r_LaBr都要加上1*cm,最好不要用 3*r_LaBr(原本的数值）
    
  /*
   G4Box*  AluminumFwallLayer  = new G4Box("AluminumFwallLayer", fAl_x/2, (height+f_delta_y+Bottom_thick)/2, aluminumdepth/2);
   G4Box*  AluminumFwallLayer  = new G4Box("AluminumFwallLayer", fAl_x/2,  fAl_y/2, aluminumdepth/2);
    
   G4Tubs* Spot_FAluminum_HPGe = new G4Tubs("Spot_FAluminum_HPGe", 0, r_HPGe+1*cm, aluminumdepth/2 + 0.002, 360.*deg, 360.*deg);
   G4Tubs* Spot_FAluminum_LaBr = new G4Tubs("Spot_FAluminum_LaBr", 0, r_LaBr+1*cm, aluminumdepth/2 + 0.002, 360.*deg, 360.*deg);
    
   G4SubtractionSolid* AluminumFwall_1hole = new G4SubtractionSolid("AluminumFwall_1hole", AluminumFwallLayer, Spot_FAluminum_HPGe, 0, Pos_hole1);
   G4SubtractionSolid* AluminumFwall = new G4SubtractionSolid("AluminumFwall", AluminumFwall_1hole, Spot_FAluminum_LaBr, 0, Pos_hole2);

   G4LogicalVolume* AluminumFwallLV = new G4LogicalVolume(AluminumFwall, AluminumMaterial, "AluminumFwallLV");
   */
    
   G4Box*  CopperFwallLayer  = new G4Box("CopperFwallLayer", fCu_x/2, fCu_y/2, copperdepth/2);
   G4Tubs* Spot_FCopper_HPGe = new G4Tubs("Spot_FCopper_HPGe", 0, r_HPGe+1*cm, copperdepth/2 + 0.002, 360.*deg, 360.*deg);
   G4Tubs* Spot_FCopper_LaBr = new G4Tubs("Spot_FCopper_LaBr", 0, r_LaBr+1*cm, copperdepth/2 + 0.002, 360.*deg, 360.*deg);


   G4SubtractionSolid* CopperFwall_1hole = new G4SubtractionSolid("CopperFwall_1hole",CopperFwallLayer, Spot_FCopper_HPGe, 0, Pos_hole1);
   G4SubtractionSolid* CopperFwall = new G4SubtractionSolid("CopperFwall", CopperFwall_1hole, Spot_FCopper_LaBr, 0, Pos_hole2);

   G4LogicalVolume* CopperFwallLV = new G4LogicalVolume(CopperFwall, CopperMaterial, "CopperFwallLV");
      CopperFwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x, front_dY, - copperdepth/2),
                 CopperFwallLV,
                 "CopperFwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

    

    G4Box*  Lead1FwallLayer  = new G4Box("Lead1FwallLayer", fPb_x/2, fPb_y/2, leaddepth/2);
    G4Tubs* Spot_FLead1_HPGe = new G4Tubs("Spot_FLead_HPGe", 0, r_HPGe+1*cm, leaddepth/2 + 0.002, 360.*deg, 360.*deg);
    G4Tubs* Spot_FLead1_LaBr = new G4Tubs("Spot_FLead_LaBr", 0, r_LaBr+1*cm, leaddepth/2 + 0.002, 360.*deg, 360.*deg);


    G4SubtractionSolid* Lead1Fwall_1hole = new G4SubtractionSolid("Lead1Fwall_1hole", Lead1FwallLayer, Spot_FLead1_HPGe, 0, Pos_hole1);
    G4SubtractionSolid* Lead1Fwall = new G4SubtractionSolid("Lead1Fwall", Lead1Fwall_1hole, Spot_FLead1_LaBr, 0, Pos_hole2);


    G4LogicalVolume* Lead1FwallLV = new G4LogicalVolume(Lead1Fwall, LeadMaterial, "Lead1FwallLV");
       Lead1FwallPV =
       new G4PVPlacement(
                  0,
                  G4ThreeVector(fW_x, front_dY, - copperdepth - leaddepth/2),
                  Lead1FwallLV,
                  "Lead1FwallPV",
                  worldLV,
                  false,
                  0,
                  fCheckOverlaps);
    


/////////////////////

   G4Box*  BP1FwallLayer  = new G4Box("BP1FwallLayer", fBP1_x/2, fBP1_y/2, BPdepth/2);
   G4Tubs* Spot_FBP1_HPGe = new G4Tubs("Spot_FBP1_HPGe", 0, r_HPGe+1*cm, BPdepth/2 + 0.002, 360.*deg, 360.*deg);
   G4Tubs* Spot_FBP1_LaBr = new G4Tubs("Spot_FBP1_LaBr", 0, r_LaBr+1*cm, BPdepth/2 + 0.002, 360.*deg, 360.*deg);


   G4SubtractionSolid* BP1Fwall_1hole = new G4SubtractionSolid("BP1Fwall_1hole",BP1FwallLayer, Spot_FBP1_HPGe, 0, Pos_hole1);
   G4SubtractionSolid* BP1Fwall = new G4SubtractionSolid("BP1Fwall", BP1Fwall_1hole, Spot_FBP1_LaBr, 0, Pos_hole2);


   G4LogicalVolume* BP1FwallLV = new G4LogicalVolume(BP1Fwall, BPMaterial, "BP1FwallLV");
      BPFwallPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x, front_dY,- copperdepth - leaddepth - BPdepth/2),
                 BP1FwallLV,
                 "BPFwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);



/////////////////////

   G4Box*  Lead2FwallLayer  = new G4Box("Lead2FwallLayer", fPb_x/2, fPb_y/2, leaddepth/2);
   G4Tubs* Spot_FLead2_HPGe = new G4Tubs("Spot_FLead_HPGe", 0, r_HPGe+1*cm, leaddepth/2 + 0.002, 360.*deg, 360.*deg);
   G4Tubs* Spot_FLead2_LaBr = new G4Tubs("Spot_FLead_LaBr", 0, r_LaBr+1*cm, leaddepth/2 + 0.002, 360.*deg, 360.*deg);


   G4SubtractionSolid* Lead2Fwall_1hole = new G4SubtractionSolid("Lead2Fwall_1hole", Lead2FwallLayer, Spot_FLead2_HPGe, 0, Pos_hole1);
   G4SubtractionSolid* Lead2Fwall = new G4SubtractionSolid("Lead2Fwall", Lead2Fwall_1hole, Spot_FLead2_LaBr, 0, Pos_hole2);


   G4LogicalVolume* Lead2FwallLV = new G4LogicalVolume(Lead2Fwall, LeadMaterial, "Lead2FwallLV");
      Lead2FwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x, front_dY, - copperdepth -leaddepth - BPdepth - leaddepth/2),
                 Lead2FwallLV,
                 "Lead2FwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


/////////////////////
   G4Box*  BP2FwallLayer  = new G4Box("BP2FwallLayer", fBP2_x/2, fBP2_y/2, BPdepth/2);
   G4Tubs* Spot_FBP2_HPGe = new G4Tubs("Spot_FBP2_HPGe", 0, r_HPGe+1*cm, BPdepth/2 + 0.002, 360.*deg, 360.*deg);
   G4Tubs* Spot_FBP2_LaBr = new G4Tubs("Spot_FBP2_LaBr", 0, r_LaBr+1*cm, BPdepth/2 + 0.002, 360.*deg, 360.*deg);


   G4SubtractionSolid* BP2Fwall_1hole = new G4SubtractionSolid("BP2Fwall_1hole",BP2FwallLayer, Spot_FBP2_HPGe, 0, Pos_hole1);
   G4SubtractionSolid* BP2Fwall = new G4SubtractionSolid("BP2Fwall", BP2Fwall_1hole, Spot_FBP2_LaBr, 0, Pos_hole2);



   G4LogicalVolume* BP2FwallLV = new G4LogicalVolume(BP2Fwall, BPMaterial, "BP2FwallLV");
      BPFwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x, front_dY,  - copperdepth - leaddepth2 - BPdepth*3/2),
                 BP2FwallLV,
                 "BPFwallPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

    
    G4Box*  Lead3FwallLayer  = new G4Box("Lead3FwallLayer", fPb_x/2, fPb_y/2, leaddepth/2);
    G4Tubs* Spot_FLead3_HPGe = new G4Tubs("Spot_FLead_HPGe", 0, r_HPGe+1*cm, leaddepth/2 + 0.002, 360.*deg, 360.*deg);
    G4Tubs* Spot_FLead3_LaBr = new G4Tubs("Spot_FLead_LaBr", 0, r_LaBr+1*cm, leaddepth/2 + 0.002, 360.*deg, 360.*deg);

    G4SubtractionSolid* Lead3Fwall_1hole = new G4SubtractionSolid("Lead3Fwall_1hole", Lead3FwallLayer, Spot_FLead3_HPGe, 0, Pos_hole1);
    G4SubtractionSolid* Lead3Fwall = new G4SubtractionSolid("Lead3Fwall", Lead3Fwall_1hole, Spot_FLead3_LaBr, 0, Pos_hole2);


    G4LogicalVolume* Lead3FwallLV = new G4LogicalVolume(Lead3Fwall, LeadMaterial, "Lead3FwallLV");
       Lead3FwallPV =
       new G4PVPlacement(
                  0,
                  G4ThreeVector(fW_x, front_dY, - copperdepth - leaddepth2 - BPdepth*2 - leaddepth/2),
                  Lead3FwallLV,
                  "Lead3FwallPV",
                  worldLV,
                  false,
                  0,
                  fCheckOverlaps);
    
//#########################TungstenSSC ###########################  // Tungsten dimension: 6" x 7" x 12.5"
  

    G4double fW_z   =  - copperdepth - leaddepth2 - BPdepth*2 - leaddepth -Wdepth/2 - leak/2;
    G4double fW_z_1 =  - copperdepth - leaddepth2 - BPdepth*2 - leaddepth -Wdepth/2 - leak;
    

   G4Box* TungstenBlock1 = new G4Box("TungstenBlock1", W_length/2, W_height/2, Wdepth_f/2);
   G4Tubs* Spot_HPGe1 = new G4Tubs("Spot_HPGe1", 0, r_HPGe2, Wdepth_f/2 + 0.004, 360.*deg, 360.*deg);
   G4Tubs* Spot_LaBr1 = new G4Tubs("Spot_LaBr1", 0, r_LaBr2, Wdepth_f/2 + 0.004, 360.*deg, 360.*deg);
   G4SubtractionSolid* TungstenONEhole1 = new G4SubtractionSolid("TungstenONEhole1", TungstenBlock1, Spot_HPGe1, 0, G4ThreeVector(-offset_Spot - fW_x, -offset_Y , 0));
   G4SubtractionSolid* TungstenSSC1 = new G4SubtractionSolid("TungstenSSC1", TungstenONEhole1, Spot_LaBr1, 0, G4ThreeVector(+offset_Spot - fW_x, -offset_Y , 0));

    
   G4Box* TungstenBlock2 = new G4Box("TungstenBlock2", W_length/2, W_height/2 , Wdepth_b/2);
   G4Tubs* Spot_HPGe2 = new G4Tubs("Spot_HPGe2", 0, r_HPGe, Wdepth_b/2 + 0.004, 360.*deg, 360.*deg);
   G4Tubs* Spot_LaBr2 = new G4Tubs("Spot_LaBr2", 0, r_LaBr, Wdepth_b/2 + 0.004, 360.*deg, 360.*deg);
   G4SubtractionSolid* TungstenONEhole2 = new G4SubtractionSolid("TungstenONEhole2", TungstenBlock2, Spot_HPGe2, 0, G4ThreeVector(-offset_Spot - fW_x, 0, 0));
   G4SubtractionSolid* TungstenSSC2 = new G4SubtractionSolid("TungstenSSC2", TungstenONEhole2, Spot_LaBr2, 0, G4ThreeVector(+offset_Spot - fW_x, 0, 0));

    
   G4UnionSolid* TungstenSSC = new G4UnionSolid("TungstenSSC", TungstenSSC1, TungstenSSC2, 0, G4ThreeVector(0, 0, 2*25.4*mm));
   G4LogicalVolume* TungstenSSCLV = new G4LogicalVolume(TungstenSSC, TungstenMaterial, "TungstenSSCLV");


      TungstenSSCPV = 
      new G4PVPlacement(
                 0,
                G4ThreeVector(fW_x, 0, fW_z_1 - Wdepth_b/2),
                //G4ThreeVector(fW_x, 0, fW_z_1),
                 TungstenSSCLV,
                 "TungstenSSCPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);
 
//#########################Virtual detector ###########################
// All Hole or Leak is related with the Virtual detector to simulate the radiation background

 G4LogicalVolume* AirHole_LaBrLV = new G4LogicalVolume(Spot_LaBr1, WorldMaterial, "AirHole_LaBrLV");

      AirHole_LaBrPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(offset_Spot, 0, fW_z_1 - Wdepth_b/2),
                 //G4ThreeVector(offset_Spot, 0, fW_z_1 - Wdepth/2 + 2*25.4*mm/2),
                 AirHole_LaBrLV,
                 "AirHole_LaBrPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


   G4LogicalVolume* AirHole_HPGeLV = new G4LogicalVolume(Spot_HPGe1, WorldMaterial, "AirHole_HPGeLV");

      AirHole_HPGePV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(-offset_Spot, 0, fW_z_1 - Wdepth_b/2),
                 //G4ThreeVector(-offset_Spot, 0, fW_z_1 - Wdepth/2 + 2*25.4*mm/2),
                 AirHole_HPGeLV,
                 "AirHole_HPGePV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

   G4Box* SSCLeak_lr = new G4Box("SSCLeak_lr", leak/2, (W_height+2*leak)/2, Wdepth/2 -0.1);
   G4LogicalVolume* SSCLeakLV_l = new G4LogicalVolume(SSCLeak_lr, WorldMaterial, "SSCLeakLV_l");
   G4LogicalVolume* SSCLeakLV_r = new G4LogicalVolume(SSCLeak_lr, WorldMaterial, "SSCLeakLV_r");


      SSCLeakPV_l =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x + (W_length+leak)/2, offset_Y , fW_z_1),
                 SSCLeakLV_l,
                 "SSCLeakPV_l",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      SSCLeakPV_r =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x - (W_length+leak)/2, offset_Y, fW_z_1),
                 SSCLeakLV_r,
                 "SSCLeakPV_r",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);




   G4Box* SSCLeak_tb = new G4Box("SSCLeak_tb", W_length/2, leak/2, Wdepth/2 -0.1);
   G4LogicalVolume* SSCLeakLV_t = new G4LogicalVolume(SSCLeak_tb, WorldMaterial, "SSCLeakLV_t");
   G4LogicalVolume* SSCLeakLV_b = new G4LogicalVolume(SSCLeak_tb, WorldMaterial, "SSCLeakLV_b");
    
   G4Box* SSCLeak_back = new G4Box("SSCLeak_back", (W_length+leak*2)/2,  (W_height+2*leak)/2, leak/2);
   G4LogicalVolume* SSCLeakLV_back = new G4LogicalVolume(SSCLeak_back, WorldMaterial, "SSCLeakLV_back");
    
    

      SSCLeakPV_t =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x, offset_Y + (W_height+leak)/2, fW_z_1 ),
                 SSCLeakLV_t,
                 "SSCLeakPV_t",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      SSCLeakPV_b =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x, offset_Y - (W_height+leak)/2, fW_z_1),
                 SSCLeakLV_b,
                 "SSCLeakPV_b",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);
    
    SSCLeakPV_back =
    new G4PVPlacement(
               0,
               G4ThreeVector(fW_x, offset_Y , fW_z_1 + Wdepth/2 + leak/2),
               SSCLeakLV_back,
               "SSCLeakPV_back",
               worldLV,
               false,
               0,
               fCheckOverlaps);

/////////////////////
    
G4double height_FLead              = 12.5*25.4*mm;
G4double support_dx                = 0.75*25.4*mm;
G4double support_Upperhole_Ylength = leak*2 + W_height;
G4double support_mid_thickness     = 1*25.4*mm;
G4double support_Lowerhole_Ylength = 2.4*25.4*mm;
G4double fPb_height                = leak*2 + W_height + support_mid_thickness*2 + support_Lowerhole_Ylength;
G4double fPb_height_dx             = (W_height/2 +leak) - fPb_height/2;
    
    
    
    
G4Box* Support_outer = new G4Box("Support_outer", (2*leak+W_length)/2, fPb_height/2, (Wdepth+leak)/2 - 0.5);
G4Box* Support_inner_upper = new G4Box("Support_inner", (2*leak + W_length)/2 + 2, (2*leak + W_height)/2 + 2, (Wdepth+leak)/2 + 0.5 );
G4Box* Support_inner_lower = new G4Box("Support_inner", (2*leak + W_length-2*support_dx)/2, support_Lowerhole_Ylength/2 , (Wdepth+leak)/2 + 0.5 );
    
G4SubtractionSolid* Support1 = new G4SubtractionSolid("Support", Support_outer, Support_inner_upper, 0, G4ThreeVector(0, -fPb_height_dx, 0));
    
G4SubtractionSolid* Support = new G4SubtractionSolid("Support", Support1, Support_inner_lower, 0, G4ThreeVector(0, -(W_height/2+leak+support_mid_thickness+support_Lowerhole_Ylength/2)-fPb_height_dx, 0));
    

  G4LogicalVolume* SupportLV = new G4LogicalVolume(Support, SteelMaterial, "SupportLV");
    
    SupportPV =
    new G4PVPlacement(
               0,
               G4ThreeVector(fW_x, fPb_height_dx ,fW_z),
               SupportLV,
               "SupportPV",
               worldLV,
               false,
               0,
               fCheckOverlaps);

    

    
    
   G4Box* FrontLeadBlock1 = new G4Box("FrontLeadBlock1", fPb_lengthL/2, fPb_height/2, (Wdepth+leak)/2 );
   G4LogicalVolume* FrontLeadBlockLV1 = new G4LogicalVolume(FrontLeadBlock1, LeadMaterial, "FrontLeadBlockLV1");

    FrontLeadBlockPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(W_length/2 + fW_x + leak + fPb_lengthL/2, fPb_height_dx, fW_z),
                 FrontLeadBlockLV1,
                 "FrontLeadBlockPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


   G4Box* FrontLeadBlock2 = new G4Box("FrontLeadBlock2", fPb_lengthR/2,  fPb_height/2, (Wdepth+leak)/2);
   G4LogicalVolume* FrontLeadBlockLV2 = new G4LogicalVolume(FrontLeadBlock2, LeadMaterial, "FrontLeadBlockLV2");
   
     FrontLeadBlockPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x - leak- (fPb_lengthR + W_length)/2, fPb_height_dx, fW_z),
                 FrontLeadBlockLV2,
                 "FrontLeadBlockPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);




   if(f_delta_y!=0)//突出来; Make the tungsten block taller and also make the concave in top shielding to reduce the radiation
   {
   G4Box* FrontLeadBlock3 = new G4Box("FrontLeadBlock3", fPb_lengthT/2, f_delta_y/2, (Wdepth+leak)/2);
   G4LogicalVolume* FrontLeadBlockLV3 = new G4LogicalVolume(FrontLeadBlock3, LeadMaterial, "FrontLeadBlockLV3");

      FrontLeadBlockPV3 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(fW_x, W_height/2 + leak + f_delta_y /2 , fW_z),
                 FrontLeadBlockLV3,
                 "FrontLeadBlockPV3",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);
  }

   
    

    
    
    

  G4double floorTheta= -atan(0.5); //这个floor指的是Top shielding 的最大那块梯形面积
  G4double floor_Zlength = 28.4*25.4*mm;
  G4double floor_Xlength = fAl_x;
  G4double floor_dX = - floor_Xlength/2 - floor_Zlength/4 + 264.12 *mm;//Unused Variable
  G4double floor_dZ = 0;//- aluminumdepth - copperdepth - leaddepth2 - 2*BPdepth;//Unused Variable
    

   G4Trap* AluminumfloorwallLayer = new G4Trap( "Aluminumfloorwall",
                                    floor_Zlength/2, floorTheta, 0, aluminumdepth/2,
                                    floor_Xlength/2, floor_Xlength/2, 0, aluminumdepth/2,
                                    floor_Xlength/2 + floor_Zlength/2, floor_Xlength/2 + floor_Zlength/2, 0);
    
//G4Trap是梯形,定义: G4Trap (const G4String &pName, G4double pDz, G4double pTheta, G4double pPhi, G4double pDy1, G4double pDx1, G4double pDx2, G4double pAlp1, G4double pDy2, G4double pDx3, G4double pDx4, G4double pAlp2)
    
   G4Trap* CopperfloorwallLayer = new G4Trap( "Copperfloorwall",
                                    floor_Zlength/2, floorTheta, 0, copperdepth/2,
                                    floor_Xlength/2, floor_Xlength/2, 0, copperdepth/2,
                                    floor_Xlength/2 + floor_Zlength/2, floor_Xlength/2 + floor_Zlength/2, 0);


   G4Trap* LeadfloorwallLayer1 = new G4Trap( "Leadfloorwall1",
                                    floor_Zlength/2, floorTheta, 0, leaddepth/2,
                                    floor_Xlength/2, floor_Xlength/2, 0, leaddepth/2,
                                    floor_Xlength/2 + floor_Zlength/2, floor_Xlength/2 + floor_Zlength/2, 0);

   G4Trap* BPfloorwallLayer = new G4Trap( "BPfloorwall",
                                    floor_Zlength/2, floorTheta, 0, BPdepth/2,
                                    floor_Xlength/2, floor_Xlength/2, 0, BPdepth/2,
                                    floor_Xlength/2 + floor_Zlength/2, floor_Xlength/2 + floor_Zlength/2, 0);


   G4Trap* LeadfloorwallLayer2 = new G4Trap( "Leadfloorwall2",
                                    floor_Zlength/2, floorTheta, 0, leaddepth/4,
                                    floor_Xlength/2, floor_Xlength/2, 0, leaddepth/4,
                                    floor_Xlength/2 + floor_Zlength/2, floor_Xlength/2 + floor_Zlength/2, 0);

//#########################Top Shielding Wall ###########################
//G4double fPb_length   = 4.0*25.4*mm; //如果不是只看Top shielding的话,记得注释掉这句话
//G4double f_delta_y    = 5*25.4*mm; //如果不是只看Top shielding的话,记得注释掉这句话

  G4double containerdepth = 0.75*25.4*mm;

  G4double T_Xlength = 33.03*25.4*mm;
  G4double T_Zlength = 36.0*25.4*mm;
  G4double T_Xlength2 = T_Xlength+T_Zlength-2*fPb_length;
 
  G4double tTheta= -atan((T_Xlength2-T_Xlength)/2/T_Zlength);

  G4double T_dX = -(T_Xlength + T_Xlength2)/4 + 4.5*25.4*mm + 264.12*mm;  //
  G4double T_dZ = - aluminumdepth - copperdepth - leaddepth - 2*BPdepth;

  G4double T_XHole   =  fPb_lengthT + 10*mm;//G4double T_XHole = 23.5*25.4*mm;
//G4double T_ZHole   =  7.6*25.4*mm;
  G4double T_ZHole   =  8.3*25.4*mm;

  G4double Roof_Leak =  0.5*25.4*mm;
  G4double Top_Leak  =  0.5*cm;

  G4double Movement=0;
  if(f_delta_y == 0) Movement = 8*25.4*mm;

    
    
   G4Trap* AluminumTwallLayer = new G4Trap( "AluminumTwallLayer",
                                    T_Zlength/2, tTheta, 0, containerdepth/2,
                                    T_Xlength/2, T_Xlength/2, 0, containerdepth/2,
                                    T_Xlength2/2, T_Xlength2/2, 0);

   // why this way define AluminumTwallLayer is wrong ?
   //G4Trap* AluminumTwallLayer = new G4Trap( "AluminumTwallLayer",
                                     //T_Zlength,
                                     //containerdepth,
                                    // T_Xlength2,
                                    // T_Xlength);
 
   G4Box* AluminumTwallHole = new G4Box("AluminumTwallHole", T_XHole/2, containerdepth/2*1.02, T_ZHole/2+0.01);

   G4SubtractionSolid* AluminumTwall = new G4SubtractionSolid("AluminumTwall", AluminumTwallLayer, AluminumTwallHole, 0, G4ThreeVector((2*T_Xlength + T_Zlength - 6*fPb_length - 2*T_XHole)/4 + 20*mm  ,0,-(T_Zlength-T_ZHole)/2-0.01 - Movement)); // Need update 20*mm


   G4LogicalVolume* AluminumTwallLV = new G4LogicalVolume(AluminumTwall, AluminumMaterial, "AluminumTwallLV");

   
     AluminumTwallPV1 =
      new G4PVPlacement(
                 0,
               //G4ThreeVector(-(T_Xlength + T_Xlength2)/4, 0,0),
                 G4ThreeVector(T_dX, offset_Y + height/2 + Top_Leak + containerdepth/2, T_Zlength/2 + T_dZ - fPb_length),
                //T_dX = -(T_Xlength + T_Xlength2)/4 + 4.5*25.4*mm + 264.12*mm;
                //T_dZ = - aluminumdepth - copperdepth - leaddepth - 2*BPdepth
                AluminumTwallLV,
                 "AluminumTwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

     AluminumTwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(T_dX, offset_Y + height/2 + Top_Leak + 3*containerdepth/2 + copperdepth+ 3*leaddepth/2 +2*BPdepth, T_Zlength/2 + T_dZ - fPb_length),
                 AluminumTwallLV,
                 "AluminumTwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


/////////////////////
   G4Trap* CopperTwallLayer = new G4Trap( "CopperTwallLayer",
                                    T_Zlength/2, tTheta, 0, copperdepth/2,
                                    T_Xlength/2, T_Xlength/2, 0, copperdepth/2,
                                    (T_Xlength+T_Zlength-2*fPb_length)/2, (T_Xlength+T_Zlength-2*fPb_length)/2, 0);


   G4Box* CopperTwallHole = new G4Box("CopperTwallHole", T_XHole/2, copperdepth/2*1.02, T_ZHole/2+0.01);

   G4SubtractionSolid* CopperTwall = new G4SubtractionSolid("CopperTwall", CopperTwallLayer, CopperTwallHole, 0, G4ThreeVector((2*T_Xlength + T_Zlength - 6*fPb_length - 2*T_XHole)/4 + 20*mm,0,-(T_Zlength-T_ZHole)/2-0.01 - Movement));


   G4LogicalVolume* CopperTwallLV = new G4LogicalVolume(CopperTwall, CopperMaterial, "CopperTwallLV");
      CopperTwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(T_dX, offset_Y + height/2 + Top_Leak + containerdepth +copperdepth/2, T_Zlength/2 + T_dZ - fPb_length),
                 CopperTwallLV,
                 "CopperTwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


/////////////////////
   G4Trap* LeadTwallLayer1 = new G4Trap( "LeadTwallLayer1",
                                    T_Zlength/2, tTheta, 0, leaddepth/2,
                                    T_Xlength/2, T_Xlength/2, 0, leaddepth/2,
                                    (T_Xlength+T_Zlength-2*fPb_length)/2, (T_Xlength+T_Zlength-2*fPb_length)/2, 0);


   G4Box* LeadTwallHole1 = new G4Box("LeadTwallHole1", T_XHole/2, leaddepth/2*1.02, T_ZHole/2 + 0.01);

   G4SubtractionSolid* LeadTwall1 = new G4SubtractionSolid("LeadTwall1", LeadTwallLayer1, LeadTwallHole1, 0, G4ThreeVector((2*T_Xlength + T_Zlength - 6*fPb_length - 2*T_XHole)/4 + 20*mm,0,-(T_Zlength-T_ZHole)/2-0.01 - Movement));


   G4LogicalVolume* LeadTwallLV1 = new G4LogicalVolume(LeadTwall1, LeadMaterial, "LeadTwallLV1");
      LeadTwallPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(T_dX, offset_Y + height/2 + Top_Leak + containerdepth + copperdepth+leaddepth/2, T_Zlength/2 + T_dZ - fPb_length),
                 LeadTwallLV1,
                 "LeadTwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);



/////////////////////
   G4Trap* LeadTwallLayer2 = new G4Trap( "LeadTwallLayer2",
                                    T_Zlength/2, tTheta, 0, leaddepth/4,
                                    T_Xlength/2, T_Xlength/2, 0, leaddepth/4,
                                    (T_Xlength+T_Zlength-2*fPb_length)/2, (T_Xlength+T_Zlength-2*fPb_length)/2, 0);


   G4Box* LeadTwallHole2 = new G4Box("LeadTwallHole2", T_XHole/2, leaddepth/4*1.02, T_ZHole/2+0.01);

   G4SubtractionSolid* LeadTwall2 = new G4SubtractionSolid("LeadTwall2", LeadTwallLayer2, LeadTwallHole2, 0, G4ThreeVector((2*T_Xlength + T_Zlength - 6*fPb_length - 2*T_XHole)/4 + 20*mm,0,-(T_Zlength-T_ZHole)/2-0.01 - Movement));


   G4LogicalVolume* LeadTwallLV2 = new G4LogicalVolume(LeadTwall2, LeadMaterial, "LeadTwallLV2");
      LeadTwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(T_dX, offset_Y + height/2 + Top_Leak + containerdepth + copperdepth + +BPdepth + 5*leaddepth/4, T_Zlength/2 + T_dZ - fPb_length),
                 LeadTwallLV2,
                 "LeadTwallPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);




/////////////////////
   G4Trap* BPTwallLayer = new G4Trap( "BPTwallLayer",
                                    T_Zlength/2, tTheta, 0, BPdepth/2,
                                    T_Xlength/2, T_Xlength/2, 0, BPdepth/2,
                                    (T_Xlength+T_Zlength-2*fPb_length)/2, (T_Xlength+T_Zlength-2*fPb_length)/2, 0);


   G4Box* BPTwallHole = new G4Box("BPTwallHole", T_XHole/2, BPdepth/2*1.02, T_ZHole/2+0.01);

   G4SubtractionSolid* BPTwall = new G4SubtractionSolid("BPTwall", BPTwallLayer, BPTwallHole, 0, G4ThreeVector((2*T_Xlength + T_Zlength - 6*fPb_length - 2*T_XHole)/4 + 20*mm,0,-(T_Zlength-T_ZHole)/2-0.01 - Movement));


   G4LogicalVolume* BPTwallLV = new G4LogicalVolume(BPTwall, BPMaterial, "BPTwallLV");

      BPTwallPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(T_dX, offset_Y + height/2 + Top_Leak + containerdepth + copperdepth+ leaddepth +BPdepth/2, T_Zlength/2 + T_dZ - fPb_length),
                 BPTwallLV,
                 "BPTwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      BPTwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(T_dX, offset_Y + height/2 + Top_Leak + containerdepth + copperdepth+ 3*leaddepth/2 +3*BPdepth/2, T_Zlength/2 + T_dZ - fPb_length),
                 BPTwallLV,
                 "BPTwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


/////////////////////
   G4Trap* AirTwallLayer = new G4Trap( "AirTwallLayer",
                                    T_Zlength/2, tTheta, 0, Top_Leak/2,
                                    T_Xlength/2, T_Xlength/2, 0, Top_Leak/2,
                                    (T_Xlength+T_Zlength-2*fPb_length)/2, (T_Xlength+T_Zlength-2*fPb_length)/2, 0);


   G4Box* AirTwallHole = new G4Box("AirTwallHole", T_XHole/2, Top_Leak/2*1.02, T_ZHole/2+0.01);

   G4SubtractionSolid* AirTwall = new G4SubtractionSolid("AirTwall", AirTwallLayer, AirTwallHole, 0, G4ThreeVector((2*T_Xlength + T_Zlength - 6*fPb_length - 2*T_XHole)/4 + 20*mm,0,-(T_Zlength-T_ZHole)/2-0.01 - Movement));


   G4LogicalVolume* AirTwallLV = new G4LogicalVolume(AirTwall,  WorldMaterial, "AirTwallLV");

      AirTwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(T_dX, offset_Y + height/2 + Top_Leak/2, T_Zlength/2 + T_dZ - fPb_length),
                 AirTwallLV,
                 "AirTwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);



//#########################Bottom Shielding Wall ###########################

    
    
    
  
    
 

  G4double bottemext_Zlength = 11.6*25.4*mm;
   
//G4double bottemext_Xlength = fPb_lengthT; // unentangle with the front shielding parameter since the front size is always changing
  G4double bottemext_Xlength = 22.78*25.4*mm;
    
   G4Box* AluminumBwallLayer_ext1 = new G4Box("AluminumBwall_ext1", bottemext_Xlength/2, aluminumdepth/2, bottemext_Zlength/2);
   G4Box* AluminumBwallLayer_ext2 = new G4Box("AluminumBwall_ext2", W_length/2 + leak , aluminumdepth*1.02/2, (Wdepth_f + Wdepth_b)/2 + 0.01);

   //G4SubtractionSolid* AluminumBwallLayer_ext = new G4SubtractionSolid("AluminumBwallLayer_ext", AluminumBwallLayer_ext1, AluminumBwallLayer_ext2, 0, G4ThreeVector(fPb_lengthR - fPb_lengthL - leak, 0, (-bottemext_Zlength + Wdepth_f + Wdepth_b)/2 - 0.01));
    G4SubtractionSolid* AluminumBwallLayer_ext = new G4SubtractionSolid("AluminumBwallLayer_ext", AluminumBwallLayer_ext1, AluminumBwallLayer_ext2, 0, G4ThreeVector(32.8*mm, 0, (-bottemext_Zlength + Wdepth_f + Wdepth_b)/2 - 0.01));
   
    
   //G4UnionSolid* AluminumBwallLayer = new G4UnionSolid("AluminumBwallLayer", AluminumfloorwallLayer, AluminumBwallLayer_ext, 0, G4ThreeVector(floor_Zlength/4, 0, -(floor_Zlength + bottemext_Zlength)/2));
    
   //G4SubtractionSolid用法: G4SubtractionSolid(”生成物名字“,减物体,”被减数物体“，旋转的参数，生成物的位置）
   //G4SubtractionSolid*是减,上一行的G4UnionSolid*是加
   //For the UnionSolid or SubtractionSolid, the first object is the one which define the position instead of the whole object define the position
    
    /*
   G4LogicalVolume* AluminumBwallLV = new G4LogicalVolume(AluminumfloorwallLayer, AluminumMaterial, "AluminumBwallLV");

      AluminumBwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(floor_dX, -height/2 - aluminumdepth/2 + offset_Y, floor_Zlength/2 + floor_dZ),
                 AluminumBwallLV,
                 "AluminumBwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);
  */


   G4Box* CopperBwallLayer_ext1 = new G4Box("CopperBwall_ext1", bottemext_Xlength/2, copperdepth/2, bottemext_Zlength/2);
   G4Box* CopperBwallLayer_ext2 = new G4Box("CopperBwall_ext2", W_length/2 + leak , copperdepth*1.02/2, (Wdepth_f + Wdepth_b)/2 + 0.01);

  // G4SubtractionSolid* CopperBwallLayer_ext = new G4SubtractionSolid("CopperBwallLayer_ext", CopperBwallLayer_ext1, CopperBwallLayer_ext2, 0, G4ThreeVector(fPb_lengthR - fPb_lengthL - leak, 0, (-bottemext_Zlength + Wdepth_f + Wdepth_b)/2 - 0.01));
    
    G4SubtractionSolid* CopperBwallLayer_ext = new G4SubtractionSolid("CopperBwallLayer_ext", CopperBwallLayer_ext1, CopperBwallLayer_ext2, 0, G4ThreeVector(32.8*mm, 0, (-bottemext_Zlength + Wdepth_f + Wdepth_b)/2 - 0.01));
    
  // G4UnionSolid* CopperBwallLayer = new G4UnionSolid("CopperBwallLayer", CopperfloorwallLayer, CopperBwallLayer_ext, 0, G4ThreeVector(floor_Zlength/4, 0, -(floor_Zlength + bottemext_Zlength)/2));


   G4LogicalVolume* CopperBwallLV = new G4LogicalVolume(CopperfloorwallLayer, CopperMaterial, "CopperBwallLV");

      CopperBwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(floor_dX, -height/2  - copperdepth/2 + offset_Y, floor_Zlength/2 + floor_dZ),
                 CopperBwallLV,
                 "CopperBwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


   G4Box* LeadBwallLayer_ext1 = new G4Box("LeadBwall_ext1", bottemext_Xlength/2, leaddepth/2, bottemext_Zlength/2);
   G4Box* LeadBwallLayer_ext2 = new G4Box("LeadBwall_ext2", W_length/2 + leak , leaddepth*1.02/2, (Wdepth_f + Wdepth_b)/2 + 0.01);

   // G4SubtractionSolid* LeadBwallLayer_ext = new G4SubtractionSolid("LeadBwallLayer_ext", LeadBwallLayer_ext1, LeadBwallLayer_ext2, 0, G4ThreeVector(fPb_lengthR - fPb_lengthL - leak, 0, (-bottemext_Zlength + Wdepth_f + Wdepth_b)/2 - 0.01));
    G4SubtractionSolid* LeadBwallLayer_ext = new G4SubtractionSolid("LeadBwallLayer_ext", LeadBwallLayer_ext1, LeadBwallLayer_ext2, 0, G4ThreeVector(32.8*mm, 0, (-bottemext_Zlength + Wdepth_f + Wdepth_b)/2 - 0.01));
    
  // G4UnionSolid* LeadBwallLayer = new G4UnionSolid("LeadBwallLayer", LeadfloorwallLayer1, LeadBwallLayer_ext, 0, G4ThreeVector(floor_Zlength/4, 0, -(floor_Zlength + bottemext_Zlength)/2));

   G4LogicalVolume* LeadBwallLV = new G4LogicalVolume(LeadfloorwallLayer1, LeadMaterial, "LeadBwallLV");

      LeadBwallPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(floor_dX, -height/2  - copperdepth - leaddepth/2 + offset_Y, floor_Zlength/2 + floor_dZ),
                 LeadBwallLV,
                 "LeadBwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      LeadBwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(floor_dX, -height/2  - copperdepth - 3*leaddepth/2 - BPdepth + offset_Y, floor_Zlength/2 + floor_dZ),
                 LeadBwallLV,
                 "LeadBwallPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


  
    

   G4Box* BPBwallLayer_ext1 = new G4Box("BPBwall_ext1", bottemext_Xlength/2, BPdepth/2, bottemext_Zlength/2);
   G4Box* BPBwallLayer_ext2 = new G4Box("BPBwall_ext2", W_length/2 + leak, BPdepth*1.02/2, (Wdepth_f + Wdepth_b)/2 + 0.01);

  // G4SubtractionSolid* BPBwallLayer_ext = new G4SubtractionSolid("BPBwallLayer_ext", BPBwallLayer_ext1, BPBwallLayer_ext2, 0, G4ThreeVector(fPb_lengthR - fPb_lengthL - leak, 0, (-bottemext_Zlength + Wdepth_f + Wdepth_b)/2 - 0.01));
    
    G4SubtractionSolid* BPBwallLayer_ext = new G4SubtractionSolid("BPBwallLayer_ext", BPBwallLayer_ext1, BPBwallLayer_ext2, 0, G4ThreeVector(32.8*mm, 0, (-bottemext_Zlength + Wdepth_f + Wdepth_b)/2 - 0.01));
    
 //  G4UnionSolid* BPBwallLayer = new G4UnionSolid("BPBwallLayer", BPfloorwallLayer, BPBwallLayer_ext, 0, G4ThreeVector(floor_Zlength/4, 0, -(floor_Zlength + bottemext_Zlength)/2));


   G4LogicalVolume* BPBwallLV = new G4LogicalVolume(BPfloorwallLayer, BPMaterial, "BPBwallLV");

      BPBwallPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(floor_dX, -height/2  - copperdepth - leaddepth -BPdepth/2 + offset_Y, floor_Zlength/2 + floor_dZ),
                 BPBwallLV,
                 "BPBwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      BPBwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(floor_dX, -height/2  - copperdepth - 2*leaddepth - 3*BPdepth/2 + offset_Y, floor_Zlength/2 + floor_dZ),
                 BPBwallLV,
                 "BPBwallPV2",
                 worldLV,
                 false,
                 0, 
                 fCheckOverlaps);


    
    
    

//#########################Right Shielding Wall ###########################

//G4double right_dX = -302.3 *mm; //302.3 *mm is 11.901575inch
  G4double right_dX = -308.4 *mm;
//G4double right_dX= -(fPb_lengthT-fPb_lengthL-leak-3.9*25.4*mm-offset_Spot);

    
  G4double rTheta=(-45)*degree;
  G4double lr_length = 25*25.4*mm; //z方向上的长度，

/*
   G4Trap* AluminumRwallLayer = new G4Trap( "AluminumRwall",
                                    lr_length/2, rTheta, 0, height/2,
                                    aluminumdepth/sqrt(2), aluminumdepth/sqrt(2), 0, height/2,
                                    aluminumdepth/sqrt(2), aluminumdepth/sqrt(2), 0);
  //above parameter need to x 2 to get the real parameter. For example, lr_length and aluminumdepth


   G4LogicalVolume* AluminumRwallLV = new G4LogicalVolume(AluminumRwallLayer, AluminumMaterial, "AluminumRwallLV");
      AluminumRwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(right_dX - lr_length/2 + (4*leaddepth + 2*copperdepth + aluminumdepth + 4*BPdepth)/sqrt(2), offset_Y, lr_length/2),
                 AluminumRwallLV,
                 "AluminumRwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

*/

   G4Trap* CopperRwallLayer = new G4Trap( "CopperRwall",
                                    lr_length/2, rTheta, 0, height/2,
                                    copperdepth/sqrt(2), copperdepth/sqrt(2), 0, height/2,
                                    copperdepth/sqrt(2), copperdepth/sqrt(2), 0);


   G4LogicalVolume* CopperRwallLV = new G4LogicalVolume(CopperRwallLayer, CopperMaterial, "CopperRwallLV");
      CopperRwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(right_dX - lr_length/2 + (4*leaddepth + copperdepth + 4*BPdepth)/sqrt(2), offset_Y, lr_length/2), 
                 CopperRwallLV,                              
                 "CopperRwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);  


   G4Trap* LeadRwallLayer = new G4Trap( "LeadRwall",
                                    lr_length/2, rTheta, 0, height/2,
                                    leaddepth/sqrt(2), leaddepth/sqrt(2), 0, height/2,
                                    leaddepth/sqrt(2), leaddepth/sqrt(2), 0);


   G4LogicalVolume* LeadRwallLV = new G4LogicalVolume(LeadRwallLayer, LeadMaterial, "LeadRwallLV");
      LeadRwallPV1 = 
      new G4PVPlacement(
                 0,
                 G4ThreeVector(right_dX - lr_length/2 + (3*leaddepth + 4*BPdepth)/sqrt(2), offset_Y, lr_length/2),
                 LeadRwallLV,                 
                 "LeadRwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);



      LeadRwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(right_dX - lr_length/2 + (1*leaddepth + 2*BPdepth)/sqrt(2), offset_Y, lr_length/2),
                 LeadRwallLV,
                 "LeadRwallPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


   G4Trap* BPRwallLayer = new G4Trap( "BPRwall",
                                    lr_length/2, rTheta, 0, height/2,
                                    BPdepth/sqrt(2), BPdepth/sqrt(2), 0, height/2,
                                    BPdepth/sqrt(2), BPdepth/sqrt(2), 0);


   G4LogicalVolume* BPRwallLV = new G4LogicalVolume(BPRwallLayer, BPMaterial, "BPRwallLV");
      BPRwallPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(right_dX - lr_length/2 + (2*leaddepth + 3*BPdepth)/sqrt(2), offset_Y, lr_length/2),
                 BPRwallLV,
                 "BPRwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);


      BPRwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(right_dX - lr_length/2 + BPdepth/sqrt(2), offset_Y, lr_length/2),
                 BPRwallLV,
                 "BPRwallPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);



//#########################Left Shielding Wall ###########################
  G4double left_dX = 264.12 *mm;

/*
   G4Box* AluminumLwallLayer = new G4Box("AluminumLwall", aluminumdepth/2, height/2, lr_length/2);
   G4LogicalVolume* AluminumLwallLV = new G4LogicalVolume(AluminumLwallLayer, AluminumMaterial, "AluminumLwallLV");
      AluminumLwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(left_dX - 2*BPdepth - 2*leaddepth - copperdepth - aluminumdepth/2, offset_Y, lr_length/2),
                 AluminumLwallLV,                               
                 "AluminumLwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps); 
*/

   G4Box* CopperLwallLayer = new G4Box("CopperLwall", copperdepth/2, height/2, lr_length/2);
   G4LogicalVolume* CopperLwallLV = new G4LogicalVolume(CopperLwallLayer, CopperMaterial, "CopperLwallLV");
      CopperLwallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(left_dX - 2*BPdepth - 2*leaddepth - copperdepth/2, offset_Y, lr_length/2),
                 CopperLwallLV,
                 "CopperLwallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);



   G4Box* LeadLwallLayer = new G4Box("LeadLwall", leaddepth/2, height/2, lr_length/2);
   G4LogicalVolume* LeadLwallLV = new G4LogicalVolume(LeadLwallLayer, LeadMaterial, "LeadLwallLV");

      LeadLwallPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(left_dX - BPdepth - leaddepth/2, offset_Y, lr_length/2),
                 LeadLwallLV,
                 "LeadLwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

      LeadLwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(left_dX - 2*BPdepth - 3*leaddepth/2, offset_Y, lr_length/2),
                 LeadLwallLV,                            
                 "LeadLwallPV1",            
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);  



   G4Box* BPLwallLayer = new G4Box("BPLwall", BPdepth/2, height/2, lr_length/2);
   G4LogicalVolume* BPLwallLV = new G4LogicalVolume(BPLwallLayer, BPMaterial, "BPLwallLV");
      BPLwallPV1 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(left_dX - BPdepth/2, offset_Y, lr_length/2),
                 BPLwallLV,
                 "BPLwallPV1",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);



      BPLwallPV2 =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(left_dX -3* BPdepth/2 - leaddepth, offset_Y, lr_length/2),
                 BPLwallLV,
                 "BPLwallPV2",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

/*
//#########################Borated Polyethylene ##########################Back wall of borated poly

   G4double Back_L = fAl_x + lr_length;

   G4Box* BPBackLayer = new G4Box("BPBackLayer", Back_L/2, height/2, BPdepth);
   G4LogicalVolume* BPBackLV = new G4LogicalVolume(BPBackLayer, BPMaterial, "BPBackLV");
      BPBackPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(left_dX - Back_L/2, offset_Y, lr_length + BPdepth),
                 BPBackLV,
                 "BPBackPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);

*/

//#########################Concrete Wall ###########################

   G4double ConcreteThick = 100*cm;
   G4Box* ConcreteWall = new G4Box("ConcreteWall", worldSizeXY*2.5, worldSizeXY*2.5, ConcreteThick/2);
   G4LogicalVolume* ConcreteWallLV = new G4LogicalVolume(ConcreteWall,ConcreteMaterial,"ConcreteWallLV");

      ConcreteWallPV =
      new G4PVPlacement(
                 0,
                 G4ThreeVector(-10*25.4, 0., 28.4*25.4 + ConcreteThick/2),
                 ConcreteWallLV,
                 "ConcreteWallPV",
                 worldLV,
                 false,
                 0,
                 fCheckOverlaps);




  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  fNofLayers = 1;
  // 
  // Sensitive detectors
  //
  auto LaBrSD 
    = new B4cCalorimeterSD("LaBrSD", "LaBrHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(LaBrSD);
  SetSensitiveDetector("fLaBrLV",LaBrSD);

  auto HPGeSD 
    = new B4cCalorimeterSD("HPGeSD", "HPGeHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(HPGeSD);
  SetSensitiveDetector("fHPGeLV",HPGeSD);

  // 
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
