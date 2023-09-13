// ----------------------------------------------------------------------------
// land_filter.cpp
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


#include "land_filter.hpp"
#include <math.h>
#include <assert.h>
#include <R_ext/Print.h>
#include <R_ext/Arith.h>
#include <R_ext/Error.h>


//-----------------------------------------------------------------------------
// sfLineSegment method implementations
//-----------------------------------------------------------------------------
bool sfLineSegment::intersects(sfLineSegment& rTruncated, const sfLineSegment& rLineSegment) const
{
  bool bIntersects = false;

  sfLineSegment NormLineSegment(rLineSegment);
  sfLineSegment Norm(*this);

  NormLineSegment.normalise();
  Norm.normalise();

  if ((Norm.A.X > NormLineSegment.B.X) || 
      (NormLineSegment.A.X > Norm.B.X))
  {
    bIntersects = false;
  }
  else
  {
    sfCoOrd Delta1(Norm.B - Norm.A);
    sfCoOrd Delta2(NormLineSegment.B - NormLineSegment.A);
    sfCoOrd Delta3(NormLineSegment.A - Norm.A);

    // Use Cramers rule
    double dxa    = Delta1.X;
    double dxc    = Delta2.X;
    double dxd    = Delta3.X;
    double dya    = Delta1.Y;
    double dyc    = Delta2.Y;
    double dyd    = Delta3.Y;
    double dNorm  = (dxc * dya) - (dxa * dyc);

    if (dNorm != 0)
    {
      double t = ((dxc * dyd) - (dxd * dyc)) / dNorm;
      double u = ((dxa * dyd) - (dya * dxd)) / dNorm;

      if ((t >= 0) && (t <= 1.0) && (u >= 0.0) && (u <= 1.0))
      {
        double x = Norm.A.X + t * dxa;
        double y = Norm.A.Y + t * dya;
        double x1 = NormLineSegment.A.X + u * dxc;
        double y1 = NormLineSegment.A.Y + u * dyc;

        bIntersects    = true;
        rTruncated.A   = A;
        rTruncated.B.X = x;
        rTruncated.B.Y = y;
      }
    }
  }

  return (bIntersects);
}


//-----------------------------------------------------------------------------
// sfBoundingBox method implementations
//-----------------------------------------------------------------------------
sfBoundingBox::sfBoundingBox(const sfLineSegment& rLineSegment)
  : TopLeft(),
    BottomRight()
{
  if (rLineSegment.A.X > rLineSegment.B.X)
  {
    TopLeft.X     = rLineSegment.B.X;
    BottomRight.X = rLineSegment.A.X;
  }
  else
  {
    TopLeft.X     = rLineSegment.A.X;
    BottomRight.X = rLineSegment.B.X;
  }

  if (rLineSegment.A.Y > rLineSegment.B.Y)
  {
    TopLeft.Y     = rLineSegment.B.Y;
    BottomRight.Y = rLineSegment.A.Y;
  }
  else
  {
    TopLeft.Y     = rLineSegment.A.Y;
    BottomRight.Y = rLineSegment.B.Y;
  }
}

//-----------------------------------------------------------------------------

sfBoundingBox::sfBoundingBox(const sfShapeFilePolygonHeader& rPolygonHeader)
  : TopLeft(sfCoOrd(rPolygonHeader.Xmin, rPolygonHeader.Ymin)),
    BottomRight(sfCoOrd(rPolygonHeader.Xmax, rPolygonHeader.Ymax))
{

}

//-----------------------------------------------------------------------------

bool sfBoundingBox::intersection(sfBoundingBox& rIntersection, const sfBoundingBox& rBoundingBox) const
{

  rIntersection.TopLeft     = TopLeft.maxCoOrd(rBoundingBox.TopLeft);
  rIntersection.BottomRight = BottomRight.minCoOrd(rBoundingBox.BottomRight);

  bool bIntersects = !((rIntersection.TopLeft.X > rIntersection.BottomRight.X) || 
                       (rIntersection.TopLeft.Y > rIntersection.BottomRight.Y));

  return (bIntersects);
}

//-----------------------------------------------------------------------------

bool sfBoundingBox::intersects(sfLineSegment& rLine) const
{
  bool bIntersects = false;

  rLine.normalise();

  double dX = rLine.B.X - rLine.A.X;
  double dY = rLine.B.Y - rLine.A.Y;

  if (::fabs(dX) >= ::fabs(dY))
  {
    double dYIntersect = ((TopLeft.X - rLine.A.X) * dY / dX) + rLine.A.Y;

    if ((dYIntersect >= TopLeft.Y) && (dYIntersect <= BottomRight.Y))
    {
      bIntersects = true;
    }
    else
    {
      dYIntersect = ((BottomRight.X - rLine.A.X) * dY / dX) + rLine.A.Y;

      if ((dYIntersect >= TopLeft.Y) && (dYIntersect <= BottomRight.Y))
      {
        bIntersects = true;
      }
      else if (dY != 0.0)
      {
        double dXIntersect = ((TopLeft.Y - rLine.A.Y) * dX / dY) + rLine.A.X;

        if ((dXIntersect >= TopLeft.X) && (dXIntersect <= BottomRight.X))
        {
          bIntersects = true;
        }
        else
        {
          dXIntersect = ((BottomRight.Y - rLine.A.Y) * dX / dY) + rLine.A.X;
          bIntersects = ((dXIntersect >= TopLeft.X) && (dXIntersect <= BottomRight.X));
        }
      }
    }
  }
  else
  {
    double dXIntersect = ((TopLeft.Y - rLine.A.Y) * dX / dY) + rLine.A.X;

    if ((dXIntersect >= TopLeft.X) && (dXIntersect <= BottomRight.X))
    {
      bIntersects = true;
    }
    else
    {
      dXIntersect = ((BottomRight.Y - rLine.A.Y) * dX / dY) + rLine.A.X;

      if ((dXIntersect >= TopLeft.X) && (dXIntersect <= BottomRight.X))
      {
        bIntersects = true;
      }
      else if (dX != 0.0)
      {
        double dYIntersect = ((TopLeft.X - rLine.A.X) * dY / dX) + rLine.A.Y;

        if ((dYIntersect >= TopLeft.Y) && (dYIntersect <= BottomRight.Y))
        {
          bIntersects = true;
        }
        else
        {
          dYIntersect = ((BottomRight.X - rLine.A.X) * dY / dX) + rLine.A.Y;
          bIntersects = ((dYIntersect >= TopLeft.Y) && (dYIntersect <= BottomRight.Y));
        }
      }
    }
  }

  return (bIntersects);
}


//-----------------------------------------------------------------------------
// The square size the lats and lons are gridded on (1 by 1 degrees).
// This number should be wholey divisible into 90.
//-----------------------------------------------------------------------------
const double sfSeaFilter::SquareSize = 1.0;


//-----------------------------------------------------------------------------
// sfPolygonAccessor method implementations
//-----------------------------------------------------------------------------
sfSeaFilter::sfPolygonAccessor::sfPolygonAccessor()
 : Record()
{
  ShapeFile     = 0;
  SavedPtr      = 0;
  Idx           = -1;
  BufferIdx     = -1;
  BufferNumber  = -1;
  BufferSize    = 0;
  Buffer        = 0;
}

//-----------------------------------------------------------------------------

sfSeaFilter::sfPolygonAccessor::sfPolygonAccessor(FILE* pShapeFile, const sfPolygonRecord& rRecord)
 : Record(rRecord)
{
  ShapeFile     = 0;
  SavedPtr      = 0;
  Idx           = -1;
  BufferIdx     = -1;
  BufferNumber  = -1;
  BufferSize    = 0;
  Buffer        = 0;

  open(pShapeFile, rRecord);
}

//-----------------------------------------------------------------------------

sfSeaFilter::sfPolygonAccessor::~sfPolygonAccessor()
{
  close();
}

//-----------------------------------------------------------------------------

void sfSeaFilter::sfPolygonAccessor::open(FILE* pShapeFile, const sfPolygonRecord& rRecord)
{
  close();

  Record        = rRecord;
  ShapeFile     = pShapeFile;
  SavedPtr      = ::ftell(ShapeFile);
  Idx           = -1;
  BufferIdx     = -1;
  BufferNumber  = -1;
  BufferSize    = (Record.Size < 65536) ? Record.Size : 65536;
  Buffer        = new sfCoOrd[BufferSize];
}

//-----------------------------------------------------------------------------

void sfSeaFilter::sfPolygonAccessor::close()
{
  if (Buffer != 0)
  {
    delete [] Buffer;
  }
  
  if (ShapeFile != 0)
  {
    ::fseek(ShapeFile, SavedPtr, SEEK_SET);
  }

  ShapeFile     = 0;
  SavedPtr      = 0;
  Idx           = -1;
  BufferIdx     = -1;
  BufferNumber  = -1;
  BufferSize    = 0;
  Buffer        = 0;
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::sfPolygonAccessor::read(sfCoOrd& rCoOrd, int nIdx) const
{
  bool bRead = false;

  if ((nIdx >= 0) && (nIdx < Record.Size))
  {
    int nBufferNumber;

    nBufferNumber = nIdx / BufferSize;
    BufferIdx     = nIdx % BufferSize;

    if (nBufferNumber != BufferNumber)
    {
      BufferNumber = nBufferNumber;

      ::fseek(ShapeFile, Record.FilePtr + (nBufferNumber * BufferSize * sizeof(sfCoOrd)), SEEK_SET);
      ::fread((void*)Buffer, BufferSize * sizeof(sfCoOrd), 1, ShapeFile);
    }

    Idx     = nIdx;
    rCoOrd  = Buffer[BufferIdx];
    bRead   = true;
  }

  return (bRead);
}

//-----------------------------------------------------------------------------

void sfSeaFilter::findEnvolvedPolygons(sfPolygonRecordPtrByLongLongMap& rResultMap, 
                                       double start_lon, 
                                       double start_lat, 
                                       double end_lon, 
                                       double end_lat) const
{
  rResultMap.clear();

  int           nStartLonSqr = lonSquare(start_lon);
  int           nStartLatSqr = latSquare(start_lat);
  int           nEndLonSqr   = lonSquare(end_lon);
  int           nEndLatSqr   = latSquare(end_lat);
  sfLineSegment rLine(sfCoOrd(start_lon, start_lat), sfCoOrd(end_lon, end_lat));
  int           nL1;
  int           nL2;
  int           nL3;
  int           nL4;

  if (nStartLonSqr <= nEndLonSqr)
  {
    nL1 = nStartLonSqr;
    nL2 = nEndLonSqr;
  }
  else
  {
    nL2 = nStartLonSqr;
    nL1 = nEndLonSqr;
  }

  if (nStartLatSqr <= nEndLatSqr)
  {
    nL3 = nStartLatSqr;
    nL4 = nEndLatSqr;
  }
  else
  {
    nL4 = nStartLatSqr;
    nL3 = nEndLatSqr;
  }

  for (int cn = nL1 ; cn <= nL2; cn++)
  {
    double dLon  = cn * SquareSize;

    for (int cm = nL3 ; cm <= nL4; cm++)
    {
      double dLat  = cm * SquareSize;
      int    nCode = encode(lonSquare(dLon), latSquare(dLat));

      for (sfPolygonRecordPtrByIntMultiMapIter Iter = PolygonMap.lower_bound(nCode) ; Iter != PolygonMap.upper_bound(nCode) ; ++Iter)
      {
        sfPolygonRecord*  pMapRecord = Iter->second;
        sfBoundingBox     BoundingBox(pMapRecord->Record);

        // Only add polygons whose bounding box intersects the line
        if (BoundingBox.intersects(rLine))
        {
          long long   nIndex = (long long)pMapRecord;
          rResultMap[nIndex] = pMapRecord;
        }
      }
    }
  }
}                                       

//-----------------------------------------------------------------------------
// Crossings algorithm point in polygon taken from 
//   http://alienryderflex.com/polygon/
//-----------------------------------------------------------------------------
bool sfSeaFilter::pointInPolygon(const sfPolygonAccessor& rAccessor, 
                                 double x, 
                                 double y) const
{
  int     i;
  int     nLimit    = rAccessor.size();
  bool    oddNodes  = false;
  sfCoOrd CoOrd_i;
  sfCoOrd CoOrd_j;

  CoOrd_i = rAccessor[0];
  
  for (i = 1 ; i < nLimit ; i++)
  {
    CoOrd_j = CoOrd_i;
    CoOrd_i = rAccessor[i];

    if (((y >= CoOrd_j.Y) && (y < CoOrd_i.Y)) || ((y >= CoOrd_i.Y) && (y < CoOrd_j.Y)))
    {
      double dX1 = CoOrd_j.X;
      double dX2 = CoOrd_i.X;
      double dY1 = CoOrd_j.Y;
      double dY2 = CoOrd_i.Y;

      // If absolute difference is greater that 180 then the line segment crosses
      // the date line, in which case, we need to adjust one of the points to
      // remove the discontinuity of the wrap around.
      if (fabs(dX2 - dX1) > 180)
      {
        if (dX1 < 0)
        {
          dX1 += 360;
        }
        else if (dX2 < 0)
        {
          dX2 += 360;
        }
        else
        {
          //Should never happen
          assert(false);
        }

        // If the x value of the test point is negative then we need to adjust 
        // that too.
        if (x < 0)
        {
          x += 360;
        }
      }

      if ((dX1 <= x) || (dX2 <= x))
      {
        // This is checking if the intersection is to the right of the point.
        if (dY1 >= dY2)
        {
          double dT;

          dY1 = CoOrd_i.Y;
          dY2 = CoOrd_j.Y;
          dT  = dX1;
          dX1 = dX2;
          dX2 = dT;
       }

        double dDX    = (dX2 - dX1);
        double dDY    = (dY2 - dY1);
        double dTest  = (x - dX1) * dDY - (y - dY1) * dDX;
        bool   bRight = (dTest > 0.0);

        oddNodes ^= bRight;
      }
    }
  }

  return oddNodes;
}

//-----------------------------------------------------------------------------

void sfSeaFilter::bigToLittle(int& nNum) const
{
  nNum = ((nNum >> 24) & 0x000000FF) + 
         ((nNum >>  8) & 0x0000FF00) + 
         ((nNum <<  8) & 0x00FF0000) + 
         ((nNum << 24) & 0xFF000000);
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::open(const char* pCoastlineShapeFile, 
                       sfShapeFileHeader& ShapeFileHeader)
{
  bool bOpened = false;

  close();

  if (pCoastlineShapeFile != 0)
  {
    ShapeFile = ::fopen(pCoastlineShapeFile, "rb");

    if (ShapeFile != 0)
    {
      // Note that this reading of the record may seem complicate but we do it
      // this way so that the reading of the structure data is correct irrespective
      // of structural aligment.
      if ((fread((void*)&ShapeFileHeader.FileCode, 4, 1, ShapeFile)     == 1) && 
          (fread((void*)&ShapeFileHeader.Unused1, 5 * 4, 1, ShapeFile)  == 1) && 
          (fread((void*)&ShapeFileHeader.FileLength, 4, 1, ShapeFile)   == 1) && 
          (fread((void*)&ShapeFileHeader.Version, 4, 1, ShapeFile)      == 1) && 
          (fread((void*)&ShapeFileHeader.ShapeType, 4, 1, ShapeFile)    == 1) && 
          (fread((void*)&ShapeFileHeader.Xmin, 8 * 8, 1, ShapeFile)     == 1))
      {
        bigToLittle(ShapeFileHeader.FileCode);
        bigToLittle(ShapeFileHeader.FileLength);

        bOpened = (ShapeFileHeader.FileCode   == 9994) && 
                  (ShapeFileHeader.Version    == 1000) && 
                  (ShapeFileHeader.ShapeType  == PolygonShapeType);

        if (bOpened)
        {
          FirstRecordPtr = ::ftell(ShapeFile);
        }
        else
        {
          close();
        }
      }
    }
  }

  return (bOpened);
}

//-----------------------------------------------------------------------------

void sfSeaFilter::close()
{
  if (ShapeFile != 0)
  {
    ::fclose(ShapeFile);

    ShapeFile         = 0;
    FirstRecordPtr    = 0;
    CurrentRecordPtr  = 0;
  }
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::findRecord(sfShapeFileRecordHeader& rRecordHeader, int nRecordPtr)
{
  bool  bMore = false;

  // Note that this reading of the record may seem complicate but we do it
  // this way so that the reading of the structure data is correct irrespective
  // of structural aligment.
  if ((::fseek(ShapeFile, nRecordPtr, SEEK_SET) == 0)                    && 
      (fread((void*)&rRecordHeader.RecordNumber, 4, 1, ShapeFile)  == 1) &&
      (fread((void*)&rRecordHeader.ContentLength, 4, 1, ShapeFile) == 1))
  {
    bigToLittle(rRecordHeader.RecordNumber);
    bigToLittle(rRecordHeader.ContentLength);

    if (rRecordHeader.RecordNumber >= 1)
    {
      CurrentRecordPtr  = nRecordPtr;
      bMore             = true;
    }
  }

  return (bMore);
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::firstRecord(sfShapeFileRecordHeader& rRecordHeader)
{
  bool bMore = findRecord(rRecordHeader, FirstRecordPtr);

  return (bMore);
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::nextRecord(sfShapeFileRecordHeader& rRecordHeader)
{
  int   nNextRecordPtr  = CurrentRecordPtr + 8 + rRecordHeader.ContentLength * 2;
  bool  bMore           = findRecord(rRecordHeader, nNextRecordPtr);

  return (bMore);
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::firstPart(sfShapeFilePolygonHeader& rPolygonHeader,
                            long& nFilePos,
                            int& nPart, 
                            int& nNumberOfPoints, 
                            int*& pPartsArray, 
                            sfPolygonAccessor* pAccessor,
                            bool bHeaderOnly)
{
  bool  bMore = false;

  nFilePos = ::ftell(ShapeFile);

  // Note that this reading of the record may seem complicate but we do it
  // this way so that the reading of the structure data is correct irrespective
  // of structural aligment.
  if ((fread((void*)&rPolygonHeader.ShapeType, 4, 1, ShapeFile) == 1) && 
      (fread((void*)&rPolygonHeader.Xmin, 4 * 8, 1, ShapeFile)  == 1) && 
      (fread((void*)&rPolygonHeader.NumParts, 4, 1, ShapeFile)  == 1) && 
      (fread((void*)&rPolygonHeader.NumPoints, 4, 1, ShapeFile) == 1) && 
      (rPolygonHeader.ShapeType == PolygonShapeType) && 
      (rPolygonHeader.NumParts  > 0)  &&
      (rPolygonHeader.NumPoints > 0))
  {
    nPart = 0;

    if (bHeaderOnly)
    {
      // Leave file pointer where we started. This is a peek into the header
      // to allow for the optimisation of not reading a polygon that does not 
      // figure in our land testing
      ::fseek(ShapeFile, nFilePos, SEEK_SET);

      nNumberOfPoints = rPolygonHeader.NumPoints;

      bMore = true;
    }
    else
    {
      if (pPartsArray != 0)
      {
        delete pPartsArray;
      }
  
      pPartsArray = new int[rPolygonHeader.NumParts];
  
      if ((pPartsArray != 0) && 
          (fread((void*)pPartsArray, rPolygonHeader.NumParts * 4, 1, ShapeFile) == 1))
      {
        nNumberOfPoints = (rPolygonHeader.NumParts == 1) ? rPolygonHeader.NumPoints - pPartsArray[0]
                                                         : pPartsArray[1] - pPartsArray[0];
  
        nFilePos = ::ftell(ShapeFile);

        if (pAccessor != 0)
        {
          sfPolygonRecord PolygonRecord(rPolygonHeader,
                                        nFilePos,
                                        nNumberOfPoints,
                                        nPart);

          pAccessor->open(ShapeFile, PolygonRecord);
        }

        bMore = true;
      }
    }
  }

  return (bMore);
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::nextPart(const sfShapeFilePolygonHeader& rPolygonHeader, 
                           long& nFilePos,
                           int& nPart, 
                           int& nNumberOfPoints, 
                           const int* pPartsArray, 
                           sfPolygonAccessor* pAccessor)
{
  bool bMore = false;

  nPart++;

  nNumberOfPoints = 0;

  if ((pPartsArray  != 0) &&
      (nPart < rPolygonHeader.NumParts))
  {
    if (nPart == rPolygonHeader.NumParts - 1)
    {
      nNumberOfPoints = (rPolygonHeader.NumParts == 1) ? rPolygonHeader.NumPoints - pPartsArray[nPart]
                                                       : pPartsArray[nPart + 1] - pPartsArray[nPart];
    }
    else
    {
      nNumberOfPoints = rPolygonHeader.NumPoints - pPartsArray[nPart];
    }

    nFilePos = ::ftell(ShapeFile);

    if (pAccessor != 0)
    {
      sfPolygonRecord PolygonRecord(rPolygonHeader,
                                    nFilePos,
                                    nNumberOfPoints,
                                    nPart);

      pAccessor->open(ShapeFile, PolygonRecord);
    }

    bMore = (nPart < rPolygonHeader.NumParts - 1);
  }

  return (bMore);
}

//-----------------------------------------------------------------------------

void sfSeaFilter::endPartEnumeration(int*& pPartsArray)
{
  if (pPartsArray != 0)
  {
    delete pPartsArray;

    pPartsArray = 0;
  }
}

//-----------------------------------------------------------------------------

sfPolygonDirection sfSeaFilter::polygonDirection(const sfPolygonAccessor& rAccessor) const
{
  double              dSum            = 0.0;
  sfPolygonDirection  nDir            = sfPolygonDirectionDegenerate;
  double              Xmin            = rAccessor.header().Xmin;
  int                 nNumberOfPoints = rAccessor.size();

  for (int i = 0 ;  i < nNumberOfPoints - 1 ; i++)
  {
    sfCoOrd Current   = rAccessor[i];
    sfCoOrd Next      = rAccessor[(i + 1) % nNumberOfPoints];
    sfCoOrd Following = rAccessor[(i + 2) % nNumberOfPoints];

    // This is to account for polygons that cross the date line if they
    // are indeed stored this way. From my reading of ESRI documentation
    // they are always split into two polygons with nothing crossing the 
    // date line but this is just in case.
    if (Current.X < Xmin)
    {
      Current.X += 360;
    }

    if (Next.X < Xmin)
    {
      Next.X += 360;
    }

    if (Following.X < Xmin)
    {
      Following.X += 360;
    }

    dSum += (Next.X - Current.X) * (Following.Y - Current.Y) - (Next.Y - Current.Y) * (Following.X - Current.X);
  }

  if (dSum > 0.0)
  {
    nDir = sfPolygonDirectionAntiClockwise;
  }
  else if (dSum < 0.0)
  {
    nDir = sfPolygonDirectionClockwise;
  }

  return (nDir);
}

//-----------------------------------------------------------------------------

void sfSeaFilter::addPolygonReference(const sfShapeFilePolygonHeader& rRecord,
                                      long nFilterPtr,
                                      int nSize,
                                      int nPart)
{
  sfPolygonRecord* pPolygonRecord = new sfPolygonRecord(rRecord,
                                                        nFilterPtr,
                                                        nSize,
                                                        nPart);

  if (pPolygonRecord != 0)                                                        
  {
    PolygonList.push_back(pPolygonRecord);

    // Put polygon definition into PolygonMap based on extent covered squares.
    // This is an optimisation to allow us to quickly find and read the needed
    // polylines from disk without having to linearly go through the file each time.
    int nLonSqrMin = lonSquare(rRecord.Xmin);
    int nLonSqrMax = lonSquare(rRecord.Xmax);
    int nLonCount  = (nLonSqrMax - nLonSqrMin) + 1;

    int nLatSqrMin = latSquare(rRecord.Ymin);
    int nLatSqrMax = latSquare(rRecord.Ymax);
    
    int nLonSquares     = (int)(360 / SquareSize);
    int nHalfLonSquares = nLonSquares >> 1;

    // This is to account for polygons which cross the date line although 
    // according to ESRI this should never be the case and such polygons
    // should be split into two but this is just in case.
    if (nLonCount <= 0)
    {
      nLonCount += nLonSquares;
    }

    if (nLatSqrMin > nLatSqrMax)
    {
      int nTemp;

      nTemp      = nLatSqrMin;
      nLatSqrMin = nLatSqrMax;
      nLatSqrMax = nTemp;
    }

    int nLonSqr = nLonSqrMin;
    
    for (int cx = 0 ; cx < nLonCount ; cx++)
    {
      for (int nLatSqr = nLatSqrMin ; nLatSqr <= nLatSqrMax ; nLatSqr++)
      {
        int nCode = encode(nLonSqr, nLatSqr);

        PolygonMap.insert(sfIntPolygonRecordPtrPair(nCode, pPolygonRecord));
      }

      nLonSqr = (nHalfLonSquares + nLonSqr + 1) % nLonSquares - nHalfLonSquares;
    }
  }
}                                      

//-----------------------------------------------------------------------------

void sfSeaFilter::endFromStartHeadingAndDistance(double& end_lon, 
                                                 double& end_lat, 
                                                 double start_lon, 
                                                 double start_lat, 
                                                 double distance_km, 
                                                 double heading_degrees) const
{
  // Convert heading from degrees to radians
  double heading = heading_degrees * M_PI / 180.0;
  
  // Earth radius in kilometers
  double radius = 6371.0;
  
  // Convert distance to radians
  double distance = distance_km / radius;
  
  // Convert latitude and longitude to radians
  start_lat *= M_PI / 180.0;
  start_lon *= M_PI / 180.0;
  
  // Calculate destination latitude
  end_lat = ::asin(::sin(start_lat) * ::cos(distance) + ::cos(start_lat) * ::sin(distance) * ::cos(heading));
  
  // Calculate destination longitude
  end_lon = start_lon + ::atan2(::sin(heading) * ::sin(distance) * ::cos(start_lat), ::cos(distance) - ::sin(start_lat) * ::sin(end_lat));
  
  // Convert destination latitude and longitude back to degrees
  end_lat *= 180.0 / M_PI;
  end_lon *= 180.0 / M_PI;
}

//-----------------------------------------------------------------------------

void sfSeaFilter::headingAndDistanceFromStartAndEnd(double& heading,
                                                    double& distance,
                                                    double start_lon, 
                                                    double start_lat, 
                                                    double end_lon, 
                                                    double end_lat) const
{
  // Convert coordinates to radians
  start_lat *= M_PI / 180.0;
  start_lon *= M_PI / 180.0;
  end_lat   *= M_PI / 180.0;
  end_lon   *= M_PI / 180.0;

  // Radius of the Earth in kilometers
  double radius = 6371.0;
  
  // Calculate the difference in longitude
  double delta_lon = end_lon - start_lon;
  double delta_lat = end_lat - start_lat;

  // Calculate heading angle using the Vincenty formula
  double y = ::sin(delta_lon) * ::cos(end_lat);
  double x = ::cos(start_lat) * ::sin(end_lat) - ::sin(start_lat) * ::cos(end_lat) * ::cos(delta_lon);
  heading  = ::atan2(y, x);

  // Convert heading angle to degrees
  heading *= 180.0 / M_PI;

  // Adjust the heading to a value between 0 and 360 degrees
  heading = ::fmod(heading + 360.0, 360.0);

  // Calculate the distance using the Haversine formula
  double a = ::pow(::sin(delta_lat / 2.0), 2) + ::cos(start_lat) * ::cos(end_lat) * ::pow(::sin(delta_lon / 2.0), 2);
  double c = 2.0 * ::atan2(::sqrt(a), ::sqrt(1.0 - a));
  
  distance =  radius * c;
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::clipEndAtLand(double& dClippedEndLon, 
                                double& dClippedEndLat, 
                                double dLonStart, 
                                double dLatStart,
                                const sfPolygonAccessor& rAccessor,
                                sfLineSegment* pIntersectingLine) const
{
  // This line clipping is done on the basis of a line drawn with longitudes and 
  // latitudes rather than true greater circle paths which is mathematically 
  // much more complicated but the approximation should be ok because the angles
  // involved will typically always be small (animals don't travel fast) and the 
  // land defining shape file has the same approximation as well (straight lines
  // drawn in an unprojected polar co-ordinate system).
  int nSize     = rAccessor.size();
  bool bClipped = false;

  if (nSize > 1)
  {
    sfLineSegment ClippedLine;
    sfLineSegment TrackLine(sfCoOrd(dLonStart, dLatStart), sfCoOrd(dClippedEndLon, dClippedEndLat));
    sfCoOrd       A = rAccessor[0];
    sfCoOrd       B;

    for (int cn = 1 ; cn < nSize ; cn++)
    {
      B = rAccessor[cn];

      sfLineSegment LandLine(A, B);

      if (TrackLine.intersects(ClippedLine, LandLine))
      {
        dClippedEndLon = ClippedLine.B.X;
        dClippedEndLat = ClippedLine.B.Y;
        TrackLine.B.X  = ClippedLine.B.X;
        TrackLine.B.Y  = ClippedLine.B.Y;
        bClipped       = true;

        if (pIntersectingLine != 0)
        {
          *pIntersectingLine = LandLine;
        }
      }

      A = B;
    }
  }

  return (bClipped);
}                                

//-----------------------------------------------------------------------------

bool sfSeaFilter::pointOnLand(double dLon, 
                              double dLat,
                              sfPolygonAccessor& rAccessor) const
{
  bool bOnLand = false;

  normalise(dLon);

  int nCode = encode(lonSquare(dLon), latSquare(dLat));

  for (sfPolygonRecordPtrByIntMultiMapConstIter Iter = PolygonMap.lower_bound(nCode) ; Iter != PolygonMap.upper_bound(nCode) ; ++Iter)
  {
    const sfPolygonRecord* pMapRecord = Iter->second;

    rAccessor.open(ShapeFile, pMapRecord[0]);

    if (pointInPolygon(rAccessor, dLon, dLat))
    {
      bOnLand = true;
      pointInPolygon(rAccessor, dLon, dLat);
      break;
    }
  }

  return (bOnLand);
}

//-----------------------------------------------------------------------------

void sfSeaFilter::flush()
{
  for (sfPolygonRecordPtrListIter Iter = PolygonList.begin() ; Iter != PolygonList.end() ; ++Iter)
  {
    sfPolygonRecord*  pRecord = *Iter;

    if (pRecord != 0)
    {
      delete pRecord;
    }
  }

  PolygonMap.clear();
  PolygonList.clear();
}

//-----------------------------------------------------------------------------

sfSeaFilter::sfSeaFilter()
 : PolygonMap(),
   PolygonList()
{
  ShapeFile         = 0;
  FirstRecordPtr    = 0;
  CurrentRecordPtr  = 0;
}

//-----------------------------------------------------------------------------

sfSeaFilter::~sfSeaFilter()
{
  flush();

  if (ShapeFile != 0)
  {
    fclose(ShapeFile);

    ShapeFile = 0;
  }
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::initialise(const char* pCoastlineShapeFile)
{
  bool bInitialised = false;

  // set true to debug loading
  if (false)
  {
    dump(pCoastlineShapeFile);
  }

  flush();

  sfShapeFileHeader ShapeFileHeader = {0};

  if (open(pCoastlineShapeFile, ShapeFileHeader))
  {
    sfShapeFileRecordHeader RecordHeader  = {0};

    bInitialised = true;

    if (firstRecord(RecordHeader))
    {
      do
      {
        int*                      pPartsArray     = 0;
        sfCoOrd*                  pPointsArray    = 0;
        int                       nPart           = 0;
        int                       nNumberOfPoints = 0;
        long                      nFilePos        = -1;
        sfShapeFilePolygonHeader  PolygonHeader;
        sfPolygonAccessor         Accessor;

        if (firstPart(PolygonHeader, 
                      nFilePos,
                      nPart, 
                      nNumberOfPoints, 
                      pPartsArray, 
                      &Accessor, 
                      false) && (nNumberOfPoints >= 3)) // Exclude the large number of 1 point records, probably places.
        {
          // Looks like the first part is always clockwise or land and the 
          // additional parts are anticlockwise or lakes. Thus we only need 
          // to add the very first polygon.
          addPolygonReference(PolygonHeader,
                              nFilePos,
                              nNumberOfPoints,
                              nPart);

          endPartEnumeration(pPartsArray);
        }
      }
      while (nextRecord(RecordHeader));
    }
  }

  return (bInitialised);
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::dump(const char* pCoastlineShapeFile)
{
  bool              bDumped         = false;
  sfShapeFileHeader ShapeFileHeader = {0};

  if (open(pCoastlineShapeFile, ShapeFileHeader))
  {
    sfShapeFileRecordHeader RecordHeader  = {0};

    Rprintf("Shape File Header\n");
    Rprintf("-----------------\n");
    Rprintf("FileCode   %d\n", ShapeFileHeader.FileCode);
    Rprintf("FileLength %d\n", ShapeFileHeader.FileLength);
    Rprintf("Version    %d\n", ShapeFileHeader.Version);
    Rprintf("ShapeType  %d\n", ShapeFileHeader.ShapeType);
    Rprintf("Xmin       %g\n", ShapeFileHeader.Xmin);
    Rprintf("Ymin       %g\n", ShapeFileHeader.Ymin);
    Rprintf("Xmax       %g\n", ShapeFileHeader.Xmax);
    Rprintf("Ymax       %g\n", ShapeFileHeader.Ymax);
    Rprintf("Zmin       %g\n", ShapeFileHeader.Zmin);
    Rprintf("Zmax       %g\n", ShapeFileHeader.Zmax);
    Rprintf("Mmin       %g\n", ShapeFileHeader.Mmin);
    Rprintf("Mmax       %g\n\n", ShapeFileHeader.Mmax);

    if (firstRecord(RecordHeader))
    {
      do
      {
        int*                      pPartsArray     = 0;
        sfCoOrd*                  pPointsArray    = 0;
        int                       nPart           = 0;
        int                       nNumberOfPoints = 0;
        long                      nFilePos        = -1;
        sfShapeFilePolygonHeader  PolygonHeader;

        Rprintf("Shape File Record Header\n");
        Rprintf("------------------------\n");
        Rprintf("RecordNumber  %d\n", RecordHeader.RecordNumber);
        Rprintf("ContentLength %d\n\n", RecordHeader.ContentLength);

        if (firstPart(PolygonHeader,
                      nFilePos, 
                      nPart, 
                      nNumberOfPoints, 
                      pPartsArray, 
                      0, 
                      true))
        {
          // There seem to be a lot of 1 point records so ignore those as they are most likely places 
          // which we are not concerned with.
          if (nNumberOfPoints > 1)
          {
            Rprintf("Shape File Polygon Header\n");
            Rprintf("-------------------------\n");
            Rprintf("ShapeType %d\n", PolygonHeader.ShapeType);
            Rprintf("Xmin      %g\n", PolygonHeader.Xmin);
            Rprintf("Ymin      %g\n", PolygonHeader.Ymin);
            Rprintf("Xmax      %g\n", PolygonHeader.Xmax);
            Rprintf("Ymax      %g\n", PolygonHeader.Ymax);
            Rprintf("NumParts  %d\n", PolygonHeader.NumParts);
            Rprintf("NumPoints %d\n\n", PolygonHeader.NumPoints);
          }

          endPartEnumeration(pPartsArray);
        }
      }
      while (nextRecord(RecordHeader));
    }

    close();

    bDumped = true;
  }

  return (bDumped);
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::pointOnLand(double dLon, 
                              double dLat) const
{
  sfPolygonAccessor rAccessor;

  // We need to normalise the angles because out calling code will not have
  // normalised angles as our state space model needs to be acting on contiguous
  // lat and lon data without dateline wrapping.
  double normalised_lon = ::fmod(dLon + 180, 360.0) - 180;
  double normalised_lat = (::fmod(dLat + 90, 180) - 90) * (1 - (2 * int((dLat + 90) / 180) & 0x1));

  return (pointOnLand(normalised_lon, 
                      normalised_lat,
                      rAccessor));
}                            

//-----------------------------------------------------------------------------

bool sfSeaFilter::landLimit(double dLonStart, 
                            double dLatStart, 
                            double dLonEnd,
                            double dLatEnd,
                            double& dLimitedLon,
                            double& dLimitedLat,
                            bool bNormalise) const
{
  bool bLimited = false;

  double dNormLonStart = dLonStart;
  double dNormLatStart = dLatStart;
  double dNormLonEnd   = dLonEnd;
  double dNormLatEnd   = dLatEnd;

  if (bNormalise)
  {
    dNormLonStart = ::fmod(dNormLonStart + 180, 360.0) - 180;
    dNormLatStart = (::fmod(dNormLatStart + 90, 180) - 90) * (1 - (2 * int((dNormLatStart + 90) / 180) & 0x1));
    dNormLonEnd   = ::fmod(dNormLonEnd + 180, 360.0) - 180;
    dNormLatEnd   = (::fmod(dNormLatEnd + 90, 180) - 90) * (1 - (2 * int((dNormLatEnd + 90) / 180) & 0x1));
  }

  dLimitedLon = dNormLonEnd;
  dLimitedLat = dNormLatEnd;

  // Limit to ocean
  sfPolygonRecordPtrByLongLongMap rPolygonsMap;

  findEnvolvedPolygons(rPolygonsMap, 
                       dLimitedLon, 
                       dLimitedLat, 
                       dNormLonStart, 
                       dNormLatStart);

  sfPolygonAccessor   rAccessor;
  sfLineSegment       rIntersectingLine;

  for (sfPolygonRecordPtrByLongLongMapIter Iter = rPolygonsMap.begin() ; Iter != rPolygonsMap.end() ; ++Iter)
  {
    const sfPolygonRecord* pRecord = Iter->second;

    rAccessor.open(ShapeFile, pRecord[0]);

    if (clipEndAtLand(dLimitedLon, 
                      dLimitedLat, 
                      dNormLonStart, 
                      dNormLatStart,
                      rAccessor, 
                      &rIntersectingLine))
    {
      bLimited = true;
    }                           
  }

  if (bLimited)
  {
    // Add a small epsilon to lat an lon to make it in sea. To do so we use a fixed
    // size coincidental epsilon scaled by 1 / sin(theta) where theta is the angle
    // between land line and track line reduced to the domain [0,pi / 2]. 
    sfLineSegment rClippedSegment(sfCoOrd(dNormLonStart, dNormLatStart), sfCoOrd(dLimitedLon, dLimitedLat));
    double        dDeltaLon   = dNormLonStart - dLimitedLon;
    double        dDeltaLat   = dNormLatStart - dLimitedLat;
    double        dLength     = ::sqrt((dDeltaLon * dDeltaLon) + (dDeltaLat * dDeltaLat));

    if (dLength > 0.0)
    {
      if (dLength < 1.0e-10)
      {
        dLimitedLon = dNormLonStart;
        dLimitedLat = dNormLatStart;
      }
      else
      {
        double        dThetaTrack = rClippedSegment.slopeAngle();
        double        dThetaLand  = rIntersectingLine.slopeAngle();
        double        dDeltaTheta = ::fmod(::fabs(dThetaTrack - dThetaLand), M_PI / 2.0);
        double        dEpsilonVal = dLength * 1.0e-3;
        double        dEpsilon    = dEpsilonVal / (::sin(dDeltaTheta) + 1.0e-6);

        if (dEpsilon > dLength)
        {
          dLimitedLon = dNormLonStart;
          dLimitedLat = dNormLatStart;
        }
        else
        {
          double        dEpsilonLon = dEpsilon * dDeltaLon / dLength;
          double        dEpsilonLat = dEpsilon * dDeltaLat / dLength;

          dLimitedLon += dEpsilonLon;
          dLimitedLat += dEpsilonLat;
        }
      }
    }

    if ((dNormLonEnd * dLimitedLon < 0.0) && (fabs(dLimitedLon) > 90.0))
    {
      // Sign has changed implying limited result has crossed the date line
      if (dNormLonEnd >= 0.0)
      {
        dLimitedLon += 180;
      }
      else
      {
        dLimitedLon -= 180;
      }
    }
  }

  return (bLimited);
}

//-----------------------------------------------------------------------------

bool sfSeaFilter::speedAndLandLimit(double dLonStart, 
                                    double dLatStart, 
                                    double dVelocityLon,
                                    double dVelocityLat,
                                    double dTimeStepS,
                                    double dMaxSpeedMpS,
                                    double& dLimitedLon,
                                    double& dLimitedLat,
                                    double& dLimitedVelocityLon, 
                                    double& dLimitedVelocityLat) const
{
  bool bLimited = false;

  if (ISNA(dTimeStepS))
  {
    Rf_error("ERROR: dTimeStepS is NA but must be a valid number. See line %d in file %s", __LINE__, __FILE__);
  }

  if (ISNA(dMaxSpeedMpS))
  {
    Rf_error("ERROR: dMaxSpeedMpS is NA but must be a valid number. See line %d in file %s", __LINE__, __FILE__);
  }

  double dSpeed = ::sqrt(dVelocityLon * dVelocityLon + dVelocityLat * dVelocityLat);

  if (dSpeed > dMaxSpeedMpS)
  {
    double dScale = dMaxSpeedMpS / dSpeed;

    dLimitedVelocityLon = dVelocityLon * dScale;
    dLimitedVelocityLat = dVelocityLat * dScale;

    bLimited = true;
  }
  else
  {
    dLimitedVelocityLon = dVelocityLon;
    dLimitedVelocityLat = dVelocityLat;
  }

  // Earth's radius in meters
  double earth_radius = 6371000.0;

  // convert to radians
  double start_lon = M_PI * dLonStart / 180.0;
  double start_lat = M_PI * dLatStart / 180.0; 

  // Convert latitudinal and longitudinal velocities to radians per second
  double lon_velocity_rad = dLimitedVelocityLon / earth_radius;
  double lat_velocity_rad = dLimitedVelocityLat / earth_radius;

  // Calculate the new latitude and longitude
  double new_lat = start_lat + lat_velocity_rad * dTimeStepS;
  double new_lon = start_lon + lon_velocity_rad * dTimeStepS / ::cos(start_lat);

  // Convert back to degrees
  new_lat *= 180.0 / M_PI;
  new_lon *= 180.0 / M_PI;

  // Normalise angles
  double normalised_new_lon = ::fmod(new_lon + 180, 360.0) - 180;
  double normalised_new_lat = (::fmod(new_lat + 90, 180) - 90) * (1 - (2 * int((new_lat + 90) / 180) & 0x1));
  double lon_offset         = new_lon - normalised_new_lon;
  double lat_offset         = new_lat - normalised_new_lat;
  double original_lon       = normalised_new_lon;
  bool   bLandLimited       = false;

  // Limit to ocean
  if (landLimit(dLonStart, 
                dLatStart, 
                normalised_new_lon,
                normalised_new_lat,
                normalised_new_lon,
                normalised_new_lat,
                false))
  {
    if ((normalised_new_lon * original_lon < 0.0) && (fabs(normalised_new_lon) > 90.0))
    {
      // Sign has changed implying limited result has crossed the date line
      if (original_lon >= 0.0)
      {
        normalised_new_lon += 180;
      }
      else
      {
        normalised_new_lon -= 180;
      }
    }

    if (::fabs(dTimeStepS) >= 1e-15)
    {
      // Need to adjust velocity to match the change
      dLimitedVelocityLon = earth_radius * ::cos(start_lat) * (normalised_new_lon - dLonStart) * M_PI / (180.0 * dTimeStepS);
      dLimitedVelocityLat = earth_radius * (normalised_new_lat - dLatStart) * M_PI / (180.0 * dTimeStepS);
    }

    bLimited = true;
  }                

  // put values back into orignal calling domain (ie. remove wrapping effects)
  normalised_new_lon += lon_offset;
  normalised_new_lat += lat_offset;

  dLimitedLon = normalised_new_lon;
  dLimitedLat = normalised_new_lat;

  return (bLimited);
}

//-----------------------------------------------------------------------------

#ifdef __R_MODULE__

IMPL_TYPE_TAG(sfSeaFilter)

//-----------------------------------------------------------------------------

EXPORT void Destroy_handler(SEXP rInstance)
{
  sfSeaFilter* pContext = (sfSeaFilter*)R_ExternalPtrAddr(rInstance);
  
  if (pContext != 0)
  {
    delete pContext;
  }
  
  R_ClearExternalPtr(rInstance);
}

//-----------------------------------------------------------------------------

EXPORT SEXP Destroy(SEXP args)
{
  SEXP rInstance;

  args = CDR(args); rInstance = CAR(args);

  PROTECT(rInstance);

  ASSERT_TYPE_TAG(rInstance, sfSeaFilter);

  Destroy_handler(rInstance);
  
  SEXP Result = Rf_allocVector(INTSXP, 1);
  
  PROTECT(Result);
  INTEGER(Result)[0] = 0;
  UNPROTECT(2);
  
  return (Result);
}

//-----------------------------------------------------------------------------

EXPORT SEXP Initialise(SEXP args)
{
  SEXP Result;
  SEXP rInstance;
  SEXP CoastlineShapeFile;
    
  args = CDR(args); CoastlineShapeFile  = CAR(args);
  
  sfSeaFilter* pContext = 0;
  bool         bOk      = true;

  PROTECT(CoastlineShapeFile);

  if (isString(CoastlineShapeFile) == 0)
  {
    Rf_error("ERROR: CoastlineShapeFile must be of type CHARSXP. See line %d in file %s", __LINE__, __FILE__);

    bOk = false;
  }

  if (bOk)
  {
    pContext = new sfSeaFilter();

    if (pContext != 0)
    {
      const char* pCoastlineShapeFile = CHAR(STRING_ELT(CoastlineShapeFile, 0));

      if (pContext->initialise(pCoastlineShapeFile))
      {
        MAKE_R_EXTERNAL_PTR(Result, pContext, Destroy_handler, sfSeaFilter);
      }
      else
      {
        Rf_error("ERROR: sfSeaFilter.initialise() failed");
      
        Result = Rf_allocVector(INTSXP, 1);

        PROTECT(Result);
        INTEGER(Result)[0] = 0;
        UNPROTECT(1);
      }
    }
    else
    {
      Rf_error("ERROR: sfSeaFilter creation failed");
      
      Result = Rf_allocVector(INTSXP, 1);
      
      PROTECT(Result);
      INTEGER(Result)[0] = 0;
      UNPROTECT(1);
    }
  }
  
  UNPROTECT(1);

  return (Result);
}

//-----------------------------------------------------------------------------

EXPORT SEXP PointOnLand(SEXP args)
{
  SEXP rInstance;
  SEXP Lon;
  SEXP Lat;
    
  args = CDR(args); rInstance = CAR(args);
  args = CDR(args); Lon       = CAR(args);
  args = CDR(args); Lat       = CAR(args);
  
  PROTECT(rInstance);

  ASSERT_TYPE_TAG(rInstance, sfSeaFilter);
  
  sfSeaFilter*  pContext = (sfSeaFilter*)R_ExternalPtrAddr(rInstance);

  if (isReal(Lon) == 0)
  {
    Rf_error("ERROR: Lon must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  
  if (isReal(Lat) == 0)
  {
    Rf_error("ERROR: Lat must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }

  if (LENGTH(Lon) != LENGTH(Lat))
  {
    Rf_error("ERROR: Lon and Lat must be the same size. See line %d in file %s", __LINE__, __FILE__);
  }

  int  nLength = LENGTH(Lon);
  SEXP Result  = Rf_allocVector(LGLSXP, nLength);

  PROTECT(Result);
  PROTECT(Lon);
  PROTECT(Lat);

  for (int cn = 0 ; cn < nLength ; cn++)
  {
    bool bOnLand = pContext->pointOnLand(REAL(Lon)[cn], 
                                         REAL(Lat)[cn]);

    LOGICAL(Result)[cn] = bOnLand ? TRUE : FALSE;                                         
  }

  UNPROTECT(4);

  return (Result);
}                        

//-----------------------------------------------------------------------------

EXPORT SEXP LandLimit(SEXP args)
{
  SEXP Result = Rf_allocVector(REALSXP, 2);
  SEXP rInstance;
  SEXP LonStart;
  SEXP LatStart;
  SEXP LonEnd;
  SEXP LatEnd;
    
  args = CDR(args); rInstance   = CAR(args);
  args = CDR(args); LonStart    = CAR(args);
  args = CDR(args); LatStart    = CAR(args);
  args = CDR(args); LonEnd      = CAR(args);
  args = CDR(args); LatEnd      = CAR(args);
  
  PROTECT(rInstance);
  PROTECT(LonStart);
  PROTECT(LatStart);
  PROTECT(LonEnd);
  PROTECT(LatEnd);

  ASSERT_TYPE_TAG(rInstance, sfSeaFilter);
  
  sfSeaFilter* pContext = (sfSeaFilter*)R_ExternalPtrAddr(rInstance);

  if (isReal(LonStart) == 0)
  {
    Rf_error("ERROR: LonStart must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  
  if (isReal(LatStart) == 0)
  {
    Rf_error("ERROR: LatStart must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  
  if (isReal(LonEnd) == 0)
  {
    Rf_error("ERROR: LonEnd must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  
  if (isReal(LatEnd) == 0)
  {
    Rf_error("ERROR: LatEnd must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  

  double dLimitedLon          = 0.0;
  double dLimitedLat          = 0.0;

  bool bLimited = pContext->landLimit(REAL(LonStart)[0], 
                                      REAL(LatStart)[0], 
                                      REAL(LonEnd)[0], 
                                      REAL(LatEnd)[0],
                                      dLimitedLon,
                                      dLimitedLat,
                                      true);

  PROTECT(Result);

  REAL(Result)[0] = dLimitedLon;
  REAL(Result)[1] = dLimitedLat;

  UNPROTECT(6);

  return (Result);
}

//-----------------------------------------------------------------------------

EXPORT SEXP SpeedAndLandLimit(SEXP args)
{
  SEXP Result = Rf_allocVector(REALSXP, 4);
  SEXP rInstance;
  SEXP LonStart;
  SEXP LatStart;
  SEXP VelocityLon;
  SEXP VelocityLat;
  SEXP TimeStepS;
  SEXP MaxSpeedMpS;
    
  args = CDR(args); rInstance   = CAR(args);
  args = CDR(args); LonStart    = CAR(args);
  args = CDR(args); LatStart    = CAR(args);
  args = CDR(args); VelocityLon = CAR(args);
  args = CDR(args); VelocityLat = CAR(args);
  args = CDR(args); TimeStepS   = CAR(args);
  args = CDR(args); MaxSpeedMpS = CAR(args);
  
  PROTECT(rInstance);
  PROTECT(LonStart);
  PROTECT(LatStart);
  PROTECT(VelocityLon);
  PROTECT(VelocityLat);
  PROTECT(TimeStepS);
  PROTECT(MaxSpeedMpS);

  ASSERT_TYPE_TAG(rInstance, sfSeaFilter);
  
  sfSeaFilter* pContext = (sfSeaFilter*)R_ExternalPtrAddr(rInstance);

  if (isReal(LonStart) == 0)
  {
    Rf_error("ERROR: LonStart must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  
  if (isReal(LatStart) == 0)
  {
    Rf_error("ERROR: LatStart must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  
  if (isReal(VelocityLon) == 0)
  {
    Rf_error("ERROR: VelocityLon must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  
  if (isReal(VelocityLat) == 0)
  {
    Rf_error("ERROR: VelocityLat must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  
  if (isReal(TimeStepS) == 0)
  {
    Rf_error("ERROR: TimeStepS must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }
  
  if (isReal(MaxSpeedMpS) == 0)
  {
    Rf_error("ERROR: MaxSpeedMpS must be of type REALSXP. See line %d in file %s", __LINE__, __FILE__);
  }

  double dLimitedLon          = 0.0;
  double dLimitedLat          = 0.0;
  double dLimitedVelocityLon  = 0.0;
  double dLimitedVelocityLat  = 0.0;

  bool bLimited = pContext->speedAndLandLimit(REAL(LonStart)[0], 
                                              REAL(LatStart)[0], 
                                              REAL(VelocityLon)[0], 
                                              REAL(VelocityLat)[0],
                                              REAL(TimeStepS)[0],
                                              REAL(MaxSpeedMpS)[0],
                                              dLimitedLon,
                                              dLimitedLat,
                                              dLimitedVelocityLon, 
                                              dLimitedVelocityLat);

  PROTECT(Result);

  REAL(Result)[0] = dLimitedLon;
  REAL(Result)[1] = dLimitedLat;
  REAL(Result)[2] = dLimitedVelocityLon;
  REAL(Result)[3] = dLimitedVelocityLat;

  UNPROTECT(8);

  return (Result);
}

//-----------------------------------------------------------------------------

static const R_CallMethodDef extMethods[] = 
{
  {"lf.destroy", (DL_FUNC)&Destroy, 1},
  {"lf.initialise", (DL_FUNC)&Initialise, 1},
  {"lf.pointOnLand", (DL_FUNC)&PointOnLand, 3},
  {"lf.speedAndLandLimit", (DL_FUNC)&SpeedAndLandLimit, 7},
  0
};

//-----------------------------------------------------------------------------

EXPORT void R_init_land_filter(DllInfo* pInfo)
{
  R_registerRoutines(pInfo, 0, 0, 0, extMethods);
}

//-----------------------------------------------------------------------------

EXPORT void R_init_libland_filter(DllInfo* pInfo)
{
  R_registerRoutines(pInfo, 0, 0, 0, extMethods);
}

#endif  //__R_MODULE__
