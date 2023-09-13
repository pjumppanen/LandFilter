// ----------------------------------------------------------------------------
// land_filter.hpp
//
// Code used in combination with Unsecent Kalman Filter to fiter out land 
// crossings of tag tracks derrived from marine species. This code provides
// the means to implement speed limiting and limiting movement to the sea
// necessary for implementing a non-linear state space model with observation 
// data to be filtered using the Unscented Kalman Filter. We make use of 
// ESRI shape files that define land-sea boundaries to to implement the sea
// limiting portion of the code.
//
// Paavo Jumppanen
// Copyright (C) 2023, CSIRO Environment
// License GNU GPLv3
// ----------------------------------------------------------------------------


#ifndef __LAND_FILTER_HPP__
#define __LAND_FILTER_HPP__


#ifdef _MSC_VER

  #ifndef _CRT_SECURE_NO_WARNINGS
    #define _CRT_SECURE_NO_WARNINGS
  #endif

  // This is here to account for VC++ percularities
  #define isnan(x) _isnan(x)
  #define EXPORT extern "C" __declspec(dllexport)

#else

  #define EXPORT extern "C" 

#endif


#include <map>
#include <list>
#include <stdio.h>
#include <math.h>


#ifdef __R_MODULE__

#define USE_RINTERNALS

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define DECL_TYPE_TAG(OBJ) \
extern const SEXP OBJ##_type_tag;

//  ----------------------------------------------------------------------------

#define IMPL_TYPE_TAG(OBJ) \
const SEXP OBJ##_type_tag = Rf_install(#OBJ"TYPE_TAG");

//  ----------------------------------------------------------------------------

#define ASSERT_TYPE_TAG(s, OBJ) \
if (!((TYPEOF(s) == EXTPTRSXP) && (R_ExternalPtrTag(s) == OBJ##_type_tag)))\
{\
  Rf_error("ERROR: rInstance must be a "#OBJ"TYPE_TAG");\
}

//  ----------------------------------------------------------------------------

#define TYPE_TAG_STRING(OBJ) #OBJ"TYPE_TAG"

//  ----------------------------------------------------------------------------

#define TYPE_TAG_OBJ(OBJ) OBJ##_type_tag

//  ----------------------------------------------------------------------------

#define MAKE_R_EXTERNAL_PTR(Result, pCONTEXT, DESTROY_FN, OBJ) \
PROTECT(Result = R_MakeExternalPtr((void*)pCONTEXT, OBJ##_type_tag, Rf_mkChar("OBJ##Class")));\
R_RegisterCFinalizer(Result, DESTROY_FN);\
UNPROTECT(1)\

#endif //__R_MODULE__


//-----------------------------------------------------------------------------
// Shape file header
//-----------------------------------------------------------------------------
// Position Field         Value       Type      Byte Order
// Byte 0   File Code     9994        Integer   Big
// Byte 4   Unused        0           Integer   Big
// Byte 8   Unused        0           Integer   Big
// Byte 12  Unused        0           Integer   Big
// Byte 16  Unused        0           Integer   Big
// Byte 20  Unused        0           Integer   Big
// Byte 24  File Length   File Length Integer   Big
// Byte 28  Version       1000        Integer   Little
// Byte 32  Shape Type    Shape Type  Integer   Little
// Byte 36  Bounding Box  Xmin        Double    Little
// Byte 44  Bounding Box  Ymin        Double    Little
// Byte 52  Bounding Box  Xmax        Double    Little
// Byte 60  Bounding Box  Ymax        Double    Little
// Byte 68* Bounding Box  Zmin        Double    Little
// Byte 76* Bounding Box  Zmax        Double    Little
// Byte 84* Bounding Box  Mmin        Double    Little
// Byte 92* Bounding Box  Mmax        Double    Little
//-----------------------------------------------------------------------------
struct sfShapeFileHeader
{
  int     FileCode;
  int     Unused1;
  int     Unused2;
  int     Unused3;
  int     Unused4;
  int     Unused5;
  int     FileLength;
  int     Version;
  int     ShapeType;
  double  Xmin;
  double  Ymin;
  double  Xmax;
  double  Ymax;
  double  Zmin;
  double  Zmax;
  double  Mmin;
  double  Mmax;
};

#define SizeOfShapeFileHeader 100

#define NullShapeType         0
#define PointShapeType        1
#define PolyLineShapeType     3
#define PolygonShapeType      5
#define MultiPointShapeType   8
#define PointZShapeType       11 
#define PolyLineZShapeType    13
#define PolygonZShapeType     15
#define MultiPointZShapeType  18
#define PointMShapeType       21
#define PolyLineMShapeType    23
#define PolygonMShapeType     25
#define MultiPointMShapeType  28
#define MultiPatchShapeType   31


//-----------------------------------------------------------------------------
// Position Field         Value       Type      Byte Order
//-----------------------------------------------------------------------------
// Byte 0   Record #      Record #    Integer   Big
// Byte 4   Content Len   Content Len Integer   Big
//-----------------------------------------------------------------------------
struct sfShapeFileRecordHeader
{
  int     RecordNumber;
  int     ContentLength;
};

#define SizeOfShapeFileRecordHeader 8


//-----------------------------------------------------------------------------
// Position Field         Value       Type      Number    Byte Order
//-----------------------------------------------------------------------------
// Byte 0   Shape Type    5           Integer   1         Little
// Byte 4   Box           Box         Double    4         Little
// Byte 36  NumParts      NumParts    Integer   1         Little
// Byte 40  NumPoints     NumPoints   Integer   1         Little
// Byte 44  Parts         Parts       Integer   NumParts  Little
// Byte X   Points        Points      Point     NumPoints Little
// Note: X = 44 + 4 * NumParts
//-----------------------------------------------------------------------------
struct sfShapeFilePolygonHeader
{
  int     ShapeType;
  double  Xmin;
  double  Ymin;
  double  Xmax;
  double  Ymax;
  int     NumParts;
  int     NumPoints;

  //---------------------------------------------------------------------------

  sfShapeFilePolygonHeader()
  {
    ShapeType = 0;
    Xmin      = 0.0;
    Ymin      = 0.0;
    Xmax      = 0.0;
    Ymax      = 0.0;
    NumParts  = 0;
    NumPoints = 0;
  };

  //---------------------------------------------------------------------------

  sfShapeFilePolygonHeader(const sfShapeFilePolygonHeader& rCopy)
  {
    ShapeType = rCopy.ShapeType;
    Xmin      = rCopy.Xmin;
    Ymin      = rCopy.Ymin;
    Xmax      = rCopy.Xmax;
    Ymax      = rCopy.Ymax;
    NumParts  = rCopy.NumParts;
    NumPoints = rCopy.NumPoints;
  };

  //---------------------------------------------------------------------------
  
  sfShapeFilePolygonHeader& operator = (const sfShapeFilePolygonHeader& rCopy)
  {
    ShapeType = rCopy.ShapeType;
    Xmin      = rCopy.Xmin;
    Ymin      = rCopy.Ymin;
    Xmax      = rCopy.Xmax;
    Ymax      = rCopy.Ymax;
    NumParts  = rCopy.NumParts;
    NumPoints = rCopy.NumPoints;

    return (*this);
  };

  //---------------------------------------------------------------------------
  
  bool operator == (const sfShapeFilePolygonHeader& rRecord)
  {
    bool bIsEqual = ((ShapeType == rRecord.ShapeType) &&
                     (Xmin      == rRecord.Xmin     ) && 
                     (Ymin      == rRecord.Ymin     ) &&
                     (Xmax      == rRecord.Xmax     ) && 
                     (Ymax      == rRecord.Ymax     ) &&
                     (NumParts  == rRecord.NumParts ) &&
                     (NumPoints == rRecord.NumPoints));

    return (bIsEqual);                     
  };
};

#define SizeOfShapeFilePolygonHeader 44

//-----------------------------------------------------------------------------

struct sfCoOrd
{
  double X;
  double Y;

  sfCoOrd();
  sfCoOrd(double dX, double dY);
  sfCoOrd(const sfCoOrd& rCopy);

  sfCoOrd         maxCoOrd(const sfCoOrd& rA) const;
  sfCoOrd         minCoOrd(const sfCoOrd& rA) const;

  sfCoOrd&        operator = (const sfCoOrd& A);
  sfCoOrd         operator + (const sfCoOrd& A);
  sfCoOrd         operator - (const sfCoOrd& A);
};

//-----------------------------------------------------------------------------

inline sfCoOrd::sfCoOrd()
{
  X = 0.0;
  Y = 0.0;
}

//-----------------------------------------------------------------------------

inline sfCoOrd::sfCoOrd(double dX, double dY)
{
  X = dX;
  Y = dY;
}

//-----------------------------------------------------------------------------

inline sfCoOrd::sfCoOrd(const sfCoOrd& rCopy)
{
  X = rCopy.X;
  Y = rCopy.Y;
}

//-----------------------------------------------------------------------------

inline sfCoOrd sfCoOrd::maxCoOrd(const sfCoOrd& rA) const
{
  sfCoOrd Max((X > rA.X) ? X : rA.Y, (Y > rA.Y) ? Y : rA.Y);

  return (Max);
}

//-----------------------------------------------------------------------------

inline sfCoOrd sfCoOrd::minCoOrd(const sfCoOrd& rA) const
{
  sfCoOrd Min((X < rA.X) ? X : rA.Y, (Y < rA.Y) ? Y : rA.Y);

  return (Min);
}

//---------------------------------------------------------------------------

inline sfCoOrd& sfCoOrd::operator = (const sfCoOrd& A)
{
  X = A.X;
  Y = A.Y;

  return (*this);
};

//---------------------------------------------------------------------------

inline sfCoOrd sfCoOrd::operator + (const sfCoOrd& A)
{
  sfCoOrd Result;

  Result.X = X + A.X;
  Result.Y = Y + A.Y;

  return (Result);
};

//---------------------------------------------------------------------------

inline sfCoOrd sfCoOrd::operator - (const sfCoOrd& A)
{
  sfCoOrd Result;

  Result.X = X - A.X;
  Result.Y = Y - A.Y;

  return (Result);
};



//-----------------------------------------------------------------------------

struct sfLineSegment
{
  sfCoOrd A;
  sfCoOrd B;

  sfLineSegment();
  sfLineSegment(const sfCoOrd& rA, const sfCoOrd& rB);
  sfLineSegment(const sfLineSegment& rCopy);

  const sfCoOrd&  ptA() const;
  const sfCoOrd&  ptB() const;
  double          slopeAngle() const;

  void            ptA(const sfCoOrd& rA);
  void            ptB(const sfCoOrd& rB);

  void            normalise();

  bool            intersects(sfLineSegment& rTruncated, const sfLineSegment& rLineSegment) const;
};

//-----------------------------------------------------------------------------

inline sfLineSegment::sfLineSegment()
 : A(),
   B()
{

}

//-----------------------------------------------------------------------------

inline sfLineSegment::sfLineSegment(const sfCoOrd& rA, const sfCoOrd& rB)
 : A(rA),
   B(rB)
{
  
}

//-----------------------------------------------------------------------------

inline sfLineSegment::sfLineSegment(const sfLineSegment& rCopy)
 : A(rCopy.A),
   B(rCopy.B)
{

}

//-----------------------------------------------------------------------------

inline const sfCoOrd& sfLineSegment::ptA() const
{
  return (A);
}

//-----------------------------------------------------------------------------

inline const sfCoOrd& sfLineSegment::ptB() const
{
  return (B);
}

//-----------------------------------------------------------------------------

inline void sfLineSegment::ptA(const sfCoOrd& rA)
{
  A = rA;
}

//-----------------------------------------------------------------------------

inline void sfLineSegment::ptB(const sfCoOrd& rB)
{
  B = rB;
}

//-----------------------------------------------------------------------------

inline double sfLineSegment::slopeAngle() const
{
  double dTheta;
  double dDeltaX = (A.X - B.X);
  double dDeltaY = (A.Y - B.Y);

  if (::fabs(dDeltaX) > ::fabs(dDeltaY))
  {
    double dSlope = dDeltaY / dDeltaX;

    dTheta = ::atan(dSlope);
  }
  else
  {
    double dInvSlope = dDeltaX / dDeltaY;

    dTheta = (M_PI / 2) - ::atan(dInvSlope);
  }

  return (dTheta);
}  

//-----------------------------------------------------------------------------

inline void sfLineSegment::normalise()
{
  if (A.X > B.X)
  {
    sfCoOrd Temp(A);

    A = B;
    B = Temp;
  }
}


//-----------------------------------------------------------------------------

struct sfBoundingBox
{
  sfCoOrd   TopLeft;
  sfCoOrd   BottomRight;

  sfBoundingBox(const sfLineSegment& rLineSegment);
  sfBoundingBox(const sfShapeFilePolygonHeader& rPolygonHeader);

  bool      intersection(sfBoundingBox& rIntersection, const sfBoundingBox& rBoundingBox) const;
  bool      intersects(sfLineSegment& rLine) const;
};


//-----------------------------------------------------------------------------

enum sfPolygonDirection
{
  sfPolygonDirectionDegenerate    = 0,
  sfPolygonDirectionClockwise     = 1,
  sfPolygonDirectionAntiClockwise = 2
};

//-----------------------------------------------------------------------------
// struct sfPolygonRecord
//-----------------------------------------------------------------------------
struct sfPolygonRecord
{
  sfShapeFilePolygonHeader  Record;
  long                      FilePtr;
  int                       Size;
  int                       Part;

  //---------------------------------------------------------------------------
  
  sfPolygonRecord()
   : Record() 
  {
    FilePtr = 0;
    Size    = 0;
    Part    = 0;
  };

  //---------------------------------------------------------------------------
  
  sfPolygonRecord(const sfShapeFilePolygonHeader& rRecord,
                  long nFilePtr,
                  int nSize,
                  int nPart)
   : Record(rRecord) 
  {
    FilePtr = nFilePtr;
    Size    = nSize;
    Part    = nPart;
  };

  //---------------------------------------------------------------------------
  
  sfPolygonRecord(const sfPolygonRecord& rCopy)
   : Record(rCopy.Record)
  {
    FilePtr = rCopy.FilePtr;
    Size    = rCopy.Size;
    Part    = rCopy.Part;
  };

  //---------------------------------------------------------------------------
  
  sfPolygonRecord& operator = (const sfPolygonRecord& rCopy)
  {
    Record  = rCopy.Record;
    FilePtr = rCopy.FilePtr;
    Size    = rCopy.Size;
    Part    = rCopy.Part;

    return (*this);
  };

  //---------------------------------------------------------------------------
  
  bool operator == (const sfPolygonRecord& rRecord)
  {
    bool bIsEqual = ((Record  == rRecord.Record ) && 
                     (FilePtr == rRecord.FilePtr) && 
                     (Size    == rRecord.Size   ) &&
                     (Part    == rRecord.Part   ));

    return (bIsEqual);
  };
};


//  ----------------------------------------------------------------------------
//  Mappings of sfPolygonRecord* by Int
//  ----------------------------------------------------------------------------
typedef std::pair<const int, sfPolygonRecord*>                sfIntPolygonRecordPtrPair;
typedef std::multimap<int, sfPolygonRecord*>                  sfPolygonRecordPtrByIntMultiMap;
typedef std::multimap<int, sfPolygonRecord*>::iterator        sfPolygonRecordPtrByIntMultiMapIter;
typedef std::multimap<int, sfPolygonRecord*>::const_iterator  sfPolygonRecordPtrByIntMultiMapConstIter;

typedef std::pair<const long long, sfPolygonRecord*>                sfLongLongPolygonRecordPtrPair;
typedef std::map<long long, sfPolygonRecord*>                       sfPolygonRecordPtrByLongLongMap;
typedef std::map<long long, sfPolygonRecord*>::iterator             sfPolygonRecordPtrByLongLongMapIter;
typedef std::map<long long, sfPolygonRecord*>::const_iterator       sfPolygonRecordPtrByLongLongMapConstIter;


//  ----------------------------------------------------------------------------
//  List of sfPolygonRecord*
//  ----------------------------------------------------------------------------
typedef std::list<sfPolygonRecord*>                           sfPolygonRecordPtrList;
typedef std::list<sfPolygonRecord*>::iterator                 sfPolygonRecordPtrListIter;
typedef std::list<sfPolygonRecord*>::const_iterator           sfPolygonRecordPtrListConstIter;


//-----------------------------------------------------------------------------
// class SeaFilter
//-----------------------------------------------------------------------------
// This class checks if points are on land described by the
// supplied shape file and if so marks them as invalid, using
// the pIsInSea as a return argument in the method filter.
//-----------------------------------------------------------------------------
// Algorithm to check if a line segment intersects land polygon, check if the
// rectangle of the line segment intersects any of the bounding rectangles for
// the polygon line segments. For each match check if the lines intersect. 
// with only one match then one point lies inside and one lies outside. With 
// two matches either both are inside or both outside, and so on. With an even
// number both must be either in or out. With odd they must be different. Use 
// the point in or out algorithm to determine which is which. 
class sfSeaFilter
{
private:
  FILE*                                   ShapeFile;
  mutable sfPolygonRecordPtrByIntMultiMap PolygonMap;
  sfPolygonRecordPtrList                  PolygonList;
  long                                    FirstRecordPtr;
  long                                    CurrentRecordPtr;

  static const double                     SquareSize;

  class sfPolygonAccessor
  {
  private:
    FILE*                   ShapeFile;
    long                    SavedPtr;
    sfPolygonRecord         Record;
    mutable int             Idx;
    mutable int             BufferIdx;
    mutable int             BufferNumber;
    mutable sfCoOrd*        Buffer;
    int                     BufferSize;

  public:
    sfPolygonAccessor();
    sfPolygonAccessor(FILE* pShapeFile, const sfPolygonRecord& rRecord);
    virtual ~sfPolygonAccessor();

    void                            open(FILE* pShapeFile, const sfPolygonRecord& rRecord);
    void                            close();

    bool                            read(sfCoOrd& rCoOrd, int nIdx) const;

    sfCoOrd                         operator [] (int nIdx) const;

    const sfShapeFilePolygonHeader& header() const;
    const sfPolygonRecord&          record() const;
    int                             size() const;
  };
  
protected:
  // We encode lattitude and longitude into a single number to give fast
  // lookup of Polygons based on membership of SquareSize by SquareSize 
  // grid element. 
  void                normalise(double& dLongitude) const;
  int                 latSquare(double dLatitude) const;
  int                 lonSquare(double dLongitude) const;
  int                 encode(int nLonSquare, int nLatSquare) const;
  void                decode(int& nLonSquare, int& nLatSquare, int nCode) const;

  void                findEnvolvedPolygons(sfPolygonRecordPtrByLongLongMap& rResultMap, 
                                           double start_lon, 
                                           double start_lat, 
                                           double end_lon, 
                                           double end_lat) const;

  bool                pointInPolygon(const sfPolygonAccessor& rIter,
                                     double tx, 
                                     double ty) const;

  void                bigToLittle(int& nNum) const;

  bool                open(const char* pCoastlineShapeFile, 
                           sfShapeFileHeader& ShapeFileHeader);

  void                close();

  bool                findRecord(sfShapeFileRecordHeader& rRecordHeader, int nRecordPtr);
  bool                firstRecord(sfShapeFileRecordHeader& rRecordHeader);
  bool                nextRecord(sfShapeFileRecordHeader& rRecordHeader);

  bool                firstPart(sfShapeFilePolygonHeader& rPolygonHeader, 
                                long& nFilePos,
                                int& nPart, 
                                int& nNumberOfPoints, 
                                int*& pPartsArray, 
                                sfPolygonAccessor* pAccessor, 
                                bool bHeaderOnly);

  bool                nextPart(const sfShapeFilePolygonHeader& rPolygonHeader, 
                               long& nFilePos,
                               int& nPart, 
                               int& nNumberOfPoints, 
                               const int* pPartsArray, 
                               sfPolygonAccessor* pAccessor);

  void                endPartEnumeration(int*& pPartsArray);

  sfPolygonDirection  polygonDirection(const sfPolygonAccessor& rAccessor) const;

  void                addPolygonReference(const sfShapeFilePolygonHeader& rRecord,
                                          long nFilePtr,
                                          int nSize,
                                          int nPart);

  void                endFromStartHeadingAndDistance(double& end_lon, 
                                                     double& end_lat, 
                                                     double start_lon, 
                                                     double start_lat, 
                                                     double distance_km, 
                                                     double heading_degrees) const;

  void                headingAndDistanceFromStartAndEnd(double& heading,
                                                        double& distance,
                                                        double start_lon, 
                                                        double start_lat, 
                                                        double end_lon, 
                                                        double end_lat) const;

  bool                clipEndAtLand(double& dClippedEndLon, 
                                    double& dClippedEndLat, 
                                    double dLonStart, 
                                    double dLatStart, 
                                    const sfPolygonAccessor& rAccessor,
                                    sfLineSegment* pIntersectingLine = 0) const;

  bool                pointOnLand(double dLon, 
                                  double dLat, 
                                  sfPolygonAccessor& rAccessor) const;

  void                flush();

public:
  sfSeaFilter();
  virtual ~sfSeaFilter();

  bool                initialise(const char* pCoastlineShapeFile);

  bool                dump(const char* pCoastlineShapeFile);

  bool                pointOnLand(double dLon, 
                                  double dLat) const;

  bool                landLimit(double dLonStart, 
                                double dLatStart, 
                                double dLonEnd,
                                double dLatEnd,
                                double& dLimitedLon,
                                double& dLimitedLat,
                                bool bNormalise) const;

  bool                speedAndLandLimit(double dLonStart, 
                                        double dLatStart, 
                                        double dVelocityLon,
                                        double dVelocityLat,
                                        double dTimeStepS,
                                        double dMaxSpeedMpS,
                                        double& dLimitedLon,
                                        double& dLimitedLat,
                                        double& dLimitedVelocityLon, 
                                        double& dLimitedVelocityLat) const;
};

//-----------------------------------------------------------------------------

inline sfCoOrd sfSeaFilter::sfPolygonAccessor::operator [] (int nIdx) const
{
  sfCoOrd Result;

  read(Result, nIdx);

  return (Result);
}

//-----------------------------------------------------------------------------

inline const sfShapeFilePolygonHeader& sfSeaFilter::sfPolygonAccessor::header() const
{
  return (Record.Record);
}

//-----------------------------------------------------------------------------

inline const sfPolygonRecord& sfSeaFilter::sfPolygonAccessor::record() const
{
  return (Record);
}

//-----------------------------------------------------------------------------

inline int sfSeaFilter::sfPolygonAccessor::size() const
{
  return (Record.Size);
}

//-----------------------------------------------------------------------------

inline void sfSeaFilter::normalise(double& dLongitude) const
{
  dLongitude = ::fmod(dLongitude + 180.0, 360.0) - 180.0;
}

//-----------------------------------------------------------------------------

inline int sfSeaFilter::latSquare(double dLatitude) const
{
  int nSquare = floor(dLatitude / SquareSize);

  return (nSquare);
}

//-----------------------------------------------------------------------------

inline int sfSeaFilter::lonSquare(double dLongitude) const
{
  int nSquare = floor(dLongitude / SquareSize);

  return (nSquare);
}

//-----------------------------------------------------------------------------

inline int sfSeaFilter::encode(int nLonSquare, int nLatSquare) const
{
  int nCode = ((nLonSquare + 512) & 0x1FF) + (((nLatSquare + 512) & 0x1FF) << 10);

  return (nCode);
}

//-----------------------------------------------------------------------------

inline void sfSeaFilter::decode(int& nLonSquare, int& nLatSquare, int nCode) const
{
  nLonSquare = ((nCode & 0x1FF) - 512);
  nLatSquare = ((nCode & 0x1FF00 >> 10) - 512);
}


#ifdef __R_MODULE__

EXPORT void Destroy_handler(SEXP rInstance);

EXPORT SEXP Destroy(SEXP args);
EXPORT SEXP Initialise(SEXP args);
EXPORT SEXP PointOnLand(SEXP args);
EXPORT SEXP LandLimit(SEXP args);
EXPORT SEXP SpeedAndLandLimit(SEXP args);

#endif  //__R_MODULE__


#endif // __LAND_FILTER_HPP__

