// FERS input Validator Function sub-system
// Outputs KML GIS readable file.
// Script written by Michael Altshuler
// University of Cape Town: ALTMIC003

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
#include <iomanip> // Include the iomanip header for setprecision

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>

using namespace std;
using namespace xercesc;

// The following two functions calculate the 3dB drop from the max gain of the antenna gain pattern
/*Start*/
double sinc_antenna_gain(double theta, double alpha, double beta, double gamma)
{
    double gain = alpha * std::pow(std::sin(beta * theta) / (beta * theta), gamma);
    return gain;
}

double find_3db_drop_angle(double alpha, double beta, double gamma)
{
    const int num_points = 1000;
    const double pi = 3.14159265358979323846;
    std::vector<double> theta(num_points);
    std::vector<double> gain(num_points);

    // Calculate gain values for each angle
    for (int i = 0; i < num_points; ++i)
    {
        theta[i] = -pi + 2.0 * pi * i / (num_points - 1);
        gain[i] = sinc_antenna_gain(theta[i], alpha, beta, gamma);
    }

    // Find the maximum gain
    double max_gain = *std::max_element(gain.begin() + num_points / 2, gain.end());

    // Convert the maximum gain to dB
    double max_gain_dB = 10.0 * std::log10(max_gain);

    // Find the target gain (3dB drop)
    double target_gain_dB = max_gain_dB - 3.0;
    double target_gain = std::pow(10.0, target_gain_dB / 10.0);

    // Find the angle index where the gain is closest to the target gain
    // Considering only positive angles (from 0 to pi)
    int idx = std::distance(gain.begin() + num_points / 2, std::min_element(gain.begin() + num_points / 2, gain.end(),
                                                                            [target_gain](double a, double b)
                                                                            { return std::abs(a - target_gain) < std::abs(b - target_gain); }));

    // Get the angle at which the 3dB drop occurs
    double angle_3dB_drop = theta[idx + num_points / 2];

    // Convert the angle to degrees
    double angle_3dB_drop_deg = angle_3dB_drop * 180.0 / pi;

    return angle_3dB_drop_deg;
}
/*End*/

// Function returns coordinates from positionWayPoint values
std::string getCoordinatesFromPositionWaypoint(const DOMElement *positionWaypointElement, double referencedLatitude, double referencedLongitude, double referenceAltitude)
{
    const XMLCh *xTag = XMLString::transcode("x");
    const XMLCh *yTag = XMLString::transcode("y");
    const XMLCh *altitudeTag = XMLString::transcode("altitude");

    double x = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(xTag)->item(0)->getTextContent()));
    double y = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(yTag)->item(0)->getTextContent()));
    double altitude = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(altitudeTag)->item(0)->getTextContent()));

    double dLongitude = referencedLongitude + x / (cos(referencedLatitude * M_PI / 180) * 111319.9);
    double dLatitude = referencedLatitude + y / 111319.9;
    double dAltitudeAboveGround = altitude - referenceAltitude;

    std::stringstream coordinates;
    coordinates << std::fixed << std::setprecision(6) << dLongitude << "," << dLatitude << "," << dAltitudeAboveGround;

    return coordinates.str();
}

// Function to calculate the destination coordinate given starting coordinate, angle (in degrees), and distance (in meters)
void calculateDestinationCoordinate(double startdLatitude, double startdLongitude, double angle, double distance, double &destdLatitude, double &destdLongitude)
{
    const double R = 6371000; // Earth's radius in meters
    double d = distance / R;  // Angular distance in radians

    // Convert degrees to radians
    double startLatRad = startdLatitude * M_PI / 180;
    double startLonRad = startdLongitude * M_PI / 180;
    double angleRad = angle * M_PI / 180;

    // Calculate destination dLatitude and dLongitude in radians
    double destLatRad = asin(sin(startLatRad) * cos(d) + cos(startLatRad) * sin(d) * cos(angleRad));
    double destLonRad = startLonRad + atan2(sin(angleRad) * sin(d) * cos(startLatRad), cos(d) - sin(startLatRad) * sin(destLatRad));

    // Convert radians to degrees
    destdLatitude = destLatRad * 180 / M_PI;
    destdLongitude = destLonRad * 180 / M_PI;
}

// Function calculates the hyperbolic path and updates the dLongitude and dLatitude values. *Not a valid interpolation path at this time (04/05/2023)
void updatedLongitudedLatitudeHyperbolic(double &dLongitude, double &dLatitude, double t, double a, double b)
{
    double x_hyperbolic = a * std::cosh(t);
    double y_hyperbolic = b * std::sinh(t);

    // Assuming the start point of the hyperbolic path is at the origin
    dLongitude += x_hyperbolic / (cos(dLatitude * M_PI / 180) * 111319.9);
    dLatitude += y_hyperbolic / 111319.9;
}

// Older implementation (deprecated)
void updatedLongitudedLatitudeCubic2(double &x, double &y, double t, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    double t2 = t * t;
    double t3 = t2 * t;

    double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    double h10 = -2.0 * t3 + 3.0 * t2;
    double h01 = t3 - 2.0 * t2 + t;
    double h11 = t3 - t2;

    x = h00 * x1 + h10 * x2 + h01 * x3 + h11 * x4;
    y = h00 * y1 + h10 * y2 + h01 * y3 + h11 * y4;
}

// Function calculates the cubic path and updates the dLongitude and dLatitude values.
void updatedLongitudedLatitudeCubic(double &newdLongitude, double &newdLatitude, double t, double dLongitude1, double dLatitude1, double dLongitude4, double dLatitude4)
{
    double t2 = t * t;
    double t3 = t2 * t;

    double controlPointAngle = 45.0 * M_PI / 180.0; // 45 degrees in radians
    double controlPointDistance = 111319.9;         // Fixed distance for control points (e.g. 111319.9 meters or 1 degree)

    // Calculate control points based on fixed angle and distance
    double x2 = dLongitude1 + controlPointDistance * cos(controlPointAngle) / (cos(dLatitude1 * M_PI / 180) * 111319.9);
    double y2 = dLatitude1 + controlPointDistance * sin(controlPointAngle) / 111319.9;

    double x3 = dLongitude4 - controlPointDistance * cos(controlPointAngle) / (cos(dLatitude4 * M_PI / 180) * 111319.9);
    double y3 = dLatitude4 - controlPointDistance * sin(controlPointAngle) / 111319.9;

    // Perform cubic interpolation
    double one_minus_t = 1 - t;
    double one_minus_t2 = one_minus_t * one_minus_t;
    double one_minus_t3 = one_minus_t2 * one_minus_t;

    newdLongitude = one_minus_t3 * dLongitude1 + 3 * one_minus_t2 * t * x2 + 3 * one_minus_t * t2 * x3 + t3 * dLongitude4;
    newdLatitude = one_minus_t3 * dLatitude1 + 3 * one_minus_t2 * t * y2 + 3 * one_minus_t * t2 * y3 + t3 * dLatitude4;
}

// Function to populate antenna maps
void populateAntennaMaps(const DOMElement *element, std::map<std::string, const DOMElement *> &isotropic_antennas, std::map<std::string, const DOMElement *> &patterned_antennas, std::map<std::string, const DOMElement *> &sinc_antennas)
{
    DOMNodeList *antennaElements = element->getElementsByTagName(XMLString::transcode("antenna"));
    for (XMLSize_t i = 0; i < antennaElements->getLength(); i++)
    {
        const DOMElement *antennaElement = dynamic_cast<const DOMElement *>(antennaElements->item(i));
        const XMLCh *nameAttr = XMLString::transcode("name");
        const XMLCh *patternAttr = XMLString::transcode("pattern");
        const XMLCh *nameValue = antennaElement->getAttribute(nameAttr);
        const XMLCh *patternValue = antennaElement->getAttribute(patternAttr);
        std::string nameStr = XMLString::transcode(nameValue);

        if (XMLString::equals(patternValue, XMLString::transcode("isotropic")))
            isotropic_antennas[nameStr] = antennaElement;
        else if (XMLString::equals(patternValue, XMLString::transcode("sinc")))
            sinc_antennas[nameStr] = antennaElement;
        else
            patterned_antennas[nameStr] = antennaElement;
    }
}

// Function to check if an antenna name is isotropic
bool IsAntennaIsotropic(const std::string &antennaName, const std::map<std::string, const DOMElement *> &mapIsotropicAntennas)
{
    return mapIsotropicAntennas.find(antennaName) != mapIsotropicAntennas.end();
}

// Function to check if an antenna name is isotropic
bool IsAntennaSinc(const std::string &antennaName, const std::map<std::string, const DOMElement *> &mapPatternedAntennas)
{
    return mapPatternedAntennas.find(antennaName) != mapPatternedAntennas.end();
}

// Function that converts degrees to radians
double deg2rad(double degrees)
{
    return degrees * M_PI / 180.0;
}

// Function to generate coordinates for a circle around a given dLatitude, dLongitude and radius
std::vector<std::pair<double, double>> generate_circle_coordinates(double lat, double lon, double radius_km, int num_points = 100)
{
    std::vector<std::pair<double, double>> circle_coordinates;
    double radius_earth = 6371.0; // Earth's radius in km

    for (int i = 0; i < num_points; i++)
    {
        double bearing = deg2rad(i * 360.0 / num_points);
        double lat_rad = deg2rad(lat);
        double lon_rad = deg2rad(lon);
        double angular_distance = radius_km / radius_earth;

        double new_lat_rad = asin(sin(lat_rad) * cos(angular_distance) +
                                  cos(lat_rad) * sin(angular_distance) * cos(bearing));
        double new_lon_rad = lon_rad + atan2(sin(bearing) * sin(angular_distance) * cos(lat_rad),
                                             cos(angular_distance) - sin(lat_rad) * sin(new_lat_rad));

        double new_lat = new_lat_rad * 180.0 / M_PI;
        double new_lon = new_lon_rad * 180.0 / M_PI;

        circle_coordinates.push_back(std::make_pair(new_lat, new_lon));
    }

    return circle_coordinates;
}

/// @brief Function to to add circular antenna radion pattern to kml
/// @param kmlFile Reference to KML file
/// @param dLatitude Latitide of centre
/// @param dLongitude Longitiude of centre
/// @param dAltitude Altitude above ground
/// @param dCircleRadius Radius of plotted ring
/// @param iNumPoints Number of points to construct the ring
void AddIsotropicRadionPatternToKML(std::ofstream &kmlFile, double dLatitude, double dLongitude, double dAltitude, double dCircleRadius, int iNumPoints)
{

    std::cout << "INFO: Adding Isotropic radiation pattern to kml at long = " + std::to_string(dLongitude) + " lat = " + std::to_string(dLatitude) << std::endl;

    kmlFile << "<Placemark>\n";
    kmlFile << "    <name>Isotropic pattern range</name>\n";
    kmlFile << "    <styleUrl>#translucentPolygon</styleUrl>\n";
    std::vector<std::pair<double, double>> circle_coordinates = generate_circle_coordinates(dLatitude, dLongitude, dCircleRadius, iNumPoints);
    kmlFile << "    <Polygon>\n";
    kmlFile << "        <extrude>1</extrude>\n";
    kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
    kmlFile << "        <outerBoundaryIs>\n";
    kmlFile << "            <LinearRing>\n";
    kmlFile << "                <coordinates>\n";

    for (const auto &coord : circle_coordinates)
        kmlFile << "                    " << coord.second << "," << coord.first << "," << dAltitude << "\n";

    // Close the circle by repeating the first point
    kmlFile << "                    " << circle_coordinates[0].second << "," << circle_coordinates[0].first << "," << dAltitude << "\n";

    kmlFile << "                </coordinates>\n";
    kmlFile << "            </LinearRing>\n";
    kmlFile << "        </outerBoundaryIs>\n";
    kmlFile << "    </Polygon>\n";
    kmlFile << "</Placemark>\n";
}

// Function to return antenna element with 'sinc' pattern
const DOMElement *getAntennaElementWithSincPattern(const DOMElement *rootElement)
{
    const XMLCh *antennaTag = XMLString::transcode("antenna");
    const XMLCh *patternAttribute = XMLString::transcode("pattern");
    DOMNodeList *antennaList = rootElement->getElementsByTagName(antennaTag);

    for (XMLSize_t i = 0; i < antennaList->getLength(); ++i)
    {
        const DOMElement *currentAntennaElement = dynamic_cast<const DOMElement *>(antennaList->item(i));
        if (XMLString::equals(currentAntennaElement->getAttribute(patternAttribute), XMLString::transcode("sinc")))
        {
            return currentAntennaElement;
        }
    }
    return nullptr;
}

// Function to process the DOMElement, extract necessary data from FERSXML file, and output accordingly to KML file
bool processPlatformElement(const DOMElement *element, std::ofstream &kmlFile, double referencedLatitude, double referencedLongitude, double referenceAltitude, DOMDocument *document)
{

    // Defining constants
    const XMLCh *platformName = XMLString::transcode("name");
    const XMLCh *platformTag = XMLString::transcode("platform");
    const XMLCh *receiverTag = XMLString::transcode("receiver");
    const XMLCh *transmitterTag = XMLString::transcode("transmitter");
    const XMLCh *targetTag = XMLString::transcode("target");
    const XMLCh *positionWaypointTag = XMLString::transcode("positionwaypoint");
    const XMLCh *xTag = XMLString::transcode("x");
    const XMLCh *yTag = XMLString::transcode("y");
    const XMLCh *altitudeTag = XMLString::transcode("altitude");
    const XMLCh *motionPathTag = XMLString::transcode("motionpath");
    const XMLCh *interpolationAttr = XMLString::transcode("interpolation");
    const XMLCh *alphaTag = XMLString::transcode("alpha");
    const XMLCh *betaTag = XMLString::transcode("beta");
    const XMLCh *gammaTag = XMLString::transcode("gamma");

    // Define maps to store isotropic and patterned antennas
    std::map<std::string, const DOMElement *> isotropic_antennas;
    std::map<std::string, const DOMElement *> sinc_antennas;
    std::map<std::string, const DOMElement *> patterned_antennas;
    populateAntennaMaps(document->getDocumentElement(), isotropic_antennas, patterned_antennas, sinc_antennas);

    // First check we have a platform to work with, otherwise we dont care
    if (!XMLString::equals(element->getTagName(), platformTag))
        return true;
    // Get the positionwaypoint element
    // Check if the getElementsByTagName() function returns a valid result before using it
    DOMNodeList *positionWaypointList = element->getElementsByTagName(positionWaypointTag);
    if (positionWaypointList->getLength() == 0)
        return false;

    // Get the motionpath element where motion path = [waybpoint1 waypoint2 ...]
    const DOMElement *motionPathElement = dynamic_cast<const DOMElement *>(element->getElementsByTagName(motionPathTag)->item(0));
    // Extract the interpolation attribute
    const XMLCh *interpolation = motionPathElement->getAttribute(interpolationAttr);
    // Determine if the interpolation is linear, hyperbolic or cubic
    bool bIsLinear = (XMLString::equals(interpolation, XMLString::transcode("linear")));
    bool bIsHyperbolic = (XMLString::equals(interpolation, XMLString::transcode("hyperbolic")));
    bool bIsCubic = (XMLString::equals(interpolation, XMLString::transcode("cubic")));

    // Get a single waypoint
    const DOMElement *positionWaypointElement = dynamic_cast<const DOMElement *>(element->getElementsByTagName(positionWaypointTag)->item(0));

    // Extract the position coordinates
    double x = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(xTag)->item(0)->getTextContent()));
    double y = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(yTag)->item(0)->getTextContent()));
    double altitude = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(altitudeTag)->item(0)->getTextContent()));
    auto strPlatformName = std::string(XMLString::transcode(element->getAttribute(platformName)));
    std::cout << "INFO: Processing " + strPlatformName + " way point: x = " + std::to_string(x) + " y = " + std::to_string(y) + " alt = " + std::to_string(altitude) << std::endl;

    // Rough estimation equirectangular projection method.
    double dLongitude = referencedLongitude + x / (cos(referencedLatitude * M_PI / 180) * 111319.9);
    double dLatitude = referencedLatitude + y / 111319.9;
    double dAltitudeAboveGround = altitude - referenceAltitude;
    std::cout << "INFO: Processing " + strPlatformName + " geographic way point: long = " + std::to_string(dLongitude) + " lat = " + std::to_string(dLatitude) + " alt = " + std::to_string(dAltitudeAboveGround) << std::endl;

    // Determine the type of placemark to use
    std::string placemarkStyle;
    if (element->getElementsByTagName(receiverTag)->getLength() > 0)
        placemarkStyle = "receiver";
    else if (element->getElementsByTagName(transmitterTag)->getLength() > 0)
        placemarkStyle = "transmitter";
    else if (element->getElementsByTagName(targetTag)->getLength() > 0)
        placemarkStyle = "target";

    // Determine antenna parameters of receiver or transmitter
    bool bIsReceiver = element->getElementsByTagName(receiverTag)->getLength() > 0;
    bool bIsTransmitter = element->getElementsByTagName(transmitterTag)->getLength() > 0;
    bool bIsTarget = element->getElementsByTagName(targetTag)->getLength() > 0;

    bool bIsIsotropicPattern = false;
    bool bIsSincPattern = false;

    if (bIsReceiver)
    {
        const DOMElement *receiverElement = dynamic_cast<const DOMElement *>(element->getElementsByTagName(receiverTag)->item(0));
        const XMLCh *antennaAttr = XMLString::transcode("antenna");
        const XMLCh *antennaValue = receiverElement->getAttribute(antennaAttr);
        std::string antennaName = XMLString::transcode(antennaValue);

        bIsIsotropicPattern = IsAntennaIsotropic(antennaName, isotropic_antennas);
        bIsSincPattern = IsAntennaSinc(antennaName, isotropic_antennas);
    }
    else if (bIsTransmitter)
    {
        const DOMElement *transmitterElement = dynamic_cast<const DOMElement *>(element->getElementsByTagName(transmitterTag)->item(0));
        const XMLCh *antennaAttr = XMLString::transcode("antenna");
        const XMLCh *antennaValue = transmitterElement->getAttribute(antennaAttr);
        std::string antennaName = XMLString::transcode(antennaValue);

        bIsIsotropicPattern = IsAntennaIsotropic(antennaName, isotropic_antennas);
        bIsSincPattern = IsAntennaSinc(antennaName, isotropic_antennas);
    }
    else if (bIsTarget)
    {
        // do nothing
    }
    else
    {
        std::cout << "INFO: Processing " + strPlatformName + " does not have member receiver or transmitter";
        throw;
    }

    // If the associated pattern is isotropic, add a circular ring of radius 20 km
    if (bIsIsotropicPattern)
    {
        double dCircleRadius = 20; // Radius in km
        int iNumPoints = 100;      // Number of points to form the circle
        AddIsotropicRadionPatternToKML(kmlFile, dLatitude, dLongitude, dAltitudeAboveGround, dCircleRadius, iNumPoints);
    }
    else if (bIsSincPattern)
    {
        // Get starting oposition and orientation of antenna
        const XMLCh *startAzimuthTag = XMLString::transcode("startazimuth");
        const XMLCh *positionWaypointTag = XMLString::transcode("positionwaypoint");

        double dStartAzimuth = std::stod(XMLString::transcode(element->getElementsByTagName(startAzimuthTag)->item(0)->getTextContent()));
        const DOMElement *positionWaypointElement = dynamic_cast<DOMElement *>(element->getElementsByTagName(positionWaypointTag)->item(0));

        // Convert positionWaypoint to coordinates and parse coordinates
        std::string strCoordinates = getCoordinatesFromPositionWaypoint(positionWaypointElement, referencedLatitude, referencedLongitude, referenceAltitude);
        double dStartLatitude, dStartLongitude, dStartAltitude;

        std::istringstream coordinatesStream(strCoordinates);
        coordinatesStream >> dStartLongitude;
        coordinatesStream.ignore(1); // skip the comma
        coordinatesStream >> dStartLatitude;
        coordinatesStream.ignore(1); // skip the comma
        coordinatesStream >> dStartAltitude;

        double dArrowLength = 20000; // Adjust this value according to the desired length of the arrow
        // Adjust startAzimuth to make 0 degrees point East
        dStartAzimuth = dStartAzimuth + 180;

        // Calculate end coordinates
        double destdLatitude, destdLongitude;
        calculateDestinationCoordinate(dStartLatitude, dStartLongitude, dStartAzimuth, dArrowLength, destdLatitude, destdLongitude);

        std::stringstream endCoordinatesStream;
        endCoordinatesStream << std::fixed << std::setprecision(6) << destdLongitude << "," << destdLatitude << "," << dStartAltitude;
        std::string endCoordinates = endCoordinatesStream.str();

        // Define values for testing
        // double alpha = 1;
        // double beta = 2;
        // double gamma = 3.6;

        // Extract the antenna element with pattern="sinc"
        const DOMElement *sincAntennaElement = getAntennaElementWithSincPattern(document->getDocumentElement());

        // Extract alpha, beta, and gamma values if the antenna element was found
        double alpha = std::stod(XMLString::transcode(sincAntennaElement->getElementsByTagName(alphaTag)->item(0)->getTextContent()));
        double beta = std::stod(XMLString::transcode(sincAntennaElement->getElementsByTagName(betaTag)->item(0)->getTextContent()));
        double gamma = std::stod(XMLString::transcode(sincAntennaElement->getElementsByTagName(gammaTag)->item(0)->getTextContent()));

        double angle_3dB_drop_deg = find_3db_drop_angle(alpha, beta, gamma);

        // Calculate end coordinates for both side lines
        double sideLine1Azimuth = dStartAzimuth - angle_3dB_drop_deg;
        double sideLine2Azimuth = dStartAzimuth + angle_3dB_drop_deg;
        double sideLine1DestdLatitude, sideLine1DestdLongitude;
        double sideLine2DestdLatitude, sideLine2DestdLongitude;

        calculateDestinationCoordinate(dStartLatitude, dStartLongitude, sideLine1Azimuth, dArrowLength, sideLine1DestdLatitude, sideLine1DestdLongitude);
        calculateDestinationCoordinate(dStartLatitude, dStartLongitude, sideLine2Azimuth, dArrowLength, sideLine2DestdLatitude, sideLine2DestdLongitude);

        std::stringstream sideLine1EndCoordinatesStream, sideLine2EndCoordinatesStream;
        sideLine1EndCoordinatesStream << std::fixed << std::setprecision(6) << sideLine1DestdLongitude << "," << sideLine1DestdLatitude << "," << dStartAltitude;
        sideLine2EndCoordinatesStream << std::fixed << std::setprecision(6) << sideLine2DestdLongitude << "," << sideLine2DestdLatitude << "," << dStartAltitude;
        std::string sideLine1EndCoordinates = sideLine1EndCoordinatesStream.str();
        std::string sideLine2EndCoordinates = sideLine2EndCoordinatesStream.str();

        // Add placemarks for side lines
        for (int i = 1; i <= 2; ++i)
        {
            std::string sideLineName = "Antenna Side Line " + std::to_string(i);
            std::string sideLineEndCoordinates = (i == 1) ? sideLine1EndCoordinates : sideLine2EndCoordinates;

            kmlFile << "<Placemark>\n";
            kmlFile << "      <name>" << sideLineName << "</name>\n";
            kmlFile << "      <styleUrl>#lineStyleBlue</styleUrl>\n";
            kmlFile << "      <LineString>\n";
            kmlFile << "            <tessellate>1</tessellate>\n";
            kmlFile << "            <coordinates>\n";
            kmlFile << "            " + strCoordinates + " " + sideLineEndCoordinates + "\n";
            kmlFile << "            </coordinates>\n";
            kmlFile << "      </LineString>\n";
            kmlFile << "</Placemark>\n";
        }

        kmlFile << "<Placemark>\n";
        kmlFile << "      <name>Antenna Direction</name>\n";
        kmlFile << "      <styleUrl>#lineStyle</styleUrl>\n";
        kmlFile << "      <LineString>\n";
        kmlFile << "            <tessellate>1</tessellate>\n";
        kmlFile << "            <coordinates>\n";
        kmlFile << "            " + strCoordinates + " " + endCoordinates + "\n";
        kmlFile << "            </coordinates>\n";
        kmlFile << "      </LineString>\n";
        kmlFile << "</Placemark>\n";

        kmlFile << "<Placemark>\n";
        kmlFile << "      <name>Antenna Arrow</name>\n";
        kmlFile << "      <styleUrl>#arrowStyle</styleUrl>\n";
        kmlFile << "      <Point>\n";
        kmlFile << "          <coordinates>" + endCoordinates + "</coordinates>\n";
        kmlFile << "      </Point>\n";
        kmlFile << "      <IconStyle>\n";
        kmlFile << "          <heading>" << dStartAzimuth << "</heading>\n";
        kmlFile << "      </IconStyle>\n";
        kmlFile << "</Placemark>\n";
    }
    else if (bIsTarget)
    {
        // Do nothing
    }
    else
        std::cout << "WARNING: antenna pattern not isotropic or sinc, no pattern plotted" << std::endl;
    ;

    // Write the placemark data to the KML file
    kmlFile << "<Placemark>\n";
    kmlFile << "    <name>" << XMLString::transcode(element->getAttribute(XMLString::transcode("name"))) << "</name>\n";
    kmlFile << "    <description>" << XMLString::transcode(element->getAttribute(XMLString::transcode("description"))) << "</description>\n";

    if (bIsReceiver)
        kmlFile << "    <styleUrl>#receiver</styleUrl>\n";
    else if (bIsTransmitter)
        kmlFile << "    <styleUrl>#transmitter</styleUrl>\n";
    else if (bIsTarget)
        kmlFile << "    <styleUrl>#target</styleUrl>\n";

    // If the interpolation is linear, hyperbolic or exponential, use the gx:Track element
    if (bIsLinear || bIsHyperbolic || bIsCubic)
    {
        kmlFile << "    <gx:Track>\n";
        if (dAltitudeAboveGround > 0)
        {
            kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
            kmlFile << "        <extrude>1</extrude>\n";
        }
        else
        {
            kmlFile << "        <altitudeMode>clampToGround</altitudeMode>\n";
        }

        // Iterate through the position waypoints
        for (XMLSize_t i = 0; i < positionWaypointList->getLength(); ++i)
        {
            const DOMElement *positionWaypointElement = dynamic_cast<const DOMElement *>(positionWaypointList->item(i));

            // Extract the position coordinates
            double x = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(xTag)->item(0)->getTextContent()));
            double y = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(yTag)->item(0)->getTextContent()));
            double altitude = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(altitudeTag)->item(0)->getTextContent()));

            // Convert the position coordinates to geographic coordinates
            double dLongitude = referencedLongitude + x / (cos(referencedLatitude * M_PI / 180) * 111319.9);
            double dLatitude = referencedLatitude + y / 111319.9;
            double dAltitudeAboveGround = altitude - referenceAltitude;

            // Extract the time value
            const XMLCh *timeTag = XMLString::transcode("time");
            double time = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(timeTag)->item(0)->getTextContent()));

            // Check if interpolation is hyperbolic
            if (bIsHyperbolic)
            {
                // Calculate the hyperbolic path and update dLongitude and dLatitude values accordingly
                double a = 0.5;                                                              // Set the desired value for 'a' based on the shape of the hyperbola
                double b = 0.5;                                                              // Set the desired value for 'b' based on the shape of the hyperbola
                double t = (double)i / (positionWaypointList->getLength() - 1) * 2.0 * M_PI; // Parameter 't' varies from 0 to 2 * PI
                updatedLongitudedLatitudeHyperbolic(dLongitude, dLatitude, t, a, b);
            }

            // Check if interpolation is cubic
            if (bIsCubic && i + 1 < positionWaypointList->getLength())
            {
                // Calculate time difference between two consecutive position waypoints
                const DOMElement *nextPositionWaypointElement = dynamic_cast<const DOMElement *>(positionWaypointList->item(i + 1));
                double nextTime = std::stod(XMLString::transcode(nextPositionWaypointElement->getElementsByTagName(timeTag)->item(0)->getTextContent()));
                double time_diff = nextTime - time;

                // Extract the position coordinates for the next waypoint
                double nextX = std::stod(XMLString::transcode(nextPositionWaypointElement->getElementsByTagName(xTag)->item(0)->getTextContent()));
                double nextY = std::stod(XMLString::transcode(nextPositionWaypointElement->getElementsByTagName(yTag)->item(0)->getTextContent()));
                double nextAltitude = std::stod(XMLString::transcode(nextPositionWaypointElement->getElementsByTagName(altitudeTag)->item(0)->getTextContent()));

                // Convert the position coordinates to geographic coordinates
                double nextdLongitude = referencedLongitude + nextX / (cos(referencedLatitude * M_PI / 180) * 111319.9);
                double nextdLatitude = referencedLatitude + nextY / 111319.9;
                double nextdAltitudeAboveGround = nextAltitude - referenceAltitude;

                // Calculate control points for cubic interpolation
                double x1 = dLongitude;
                double y1 = dLatitude;
                double x4 = nextdLongitude;
                double y4 = nextdLatitude;
                double newX, newY;
                updatedLongitudedLatitudeCubic(newX, newY, 0.0, x1, y1, x4, y4); // Calculate first point on cubic curve

                int num_divisions = 100;
                for (int j = 0; j <= num_divisions; ++j)
                {
                    double t = (double)j / num_divisions;

                    double newdLongitude, newdLatitude;
                    updatedLongitudedLatitudeCubic(newdLongitude, newdLatitude, t, x1, y1, x4, y4);

                    double newdAltitudeAboveGround = dAltitudeAboveGround + t * (nextdAltitudeAboveGround - dAltitudeAboveGround);

                    kmlFile << "        <when>" << time + (double)(j * time_diff) / num_divisions << "</when>\n";
                    kmlFile << "        <gx:coord>" << newdLongitude << " " << newdLatitude << " " << newdAltitudeAboveGround << "</gx:coord>\n";
                }
            }

            else
            {
                // Write the time and coordinates to the gx:Track element
                kmlFile << "        <when>" << time << "</when>\n";
                kmlFile << "        <gx:coord>" << dLongitude << " " << dLatitude << " " << dAltitudeAboveGround << "</gx:coord>\n";
            }
        }

        kmlFile << "    </gx:Track>\n";
    }

    else
    {
        kmlFile << "    <LookAt>\n";
        kmlFile << "        <dLongitude>" << dLongitude << "</dLongitude>\n";
        kmlFile << "        <dLatitude>" << dLatitude << "</dLatitude>\n";
        kmlFile << "        <altitude>" << dAltitudeAboveGround << "</altitude>\n";
        kmlFile << "        <heading>-148.4122922628044</heading>\n";
        kmlFile << "        <tilt>40.5575073395506</tilt>\n";
        kmlFile << "        <range>500.6566641072245</range>\n";
        kmlFile << "    </LookAt>\n";

        kmlFile << "    <Point>\n";
        kmlFile << "        <coordinates>" << dLongitude << "," << dLatitude << "," << dAltitudeAboveGround << "</coordinates>\n";

        if (dAltitudeAboveGround > 0)
        {
            kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
            kmlFile << "        <extrude>1</extrude>\n";
        }
        else
            kmlFile << "        <altitudeMode>clampToGround</altitudeMode>\n";

        kmlFile << "    </Point>\n";
    }

    kmlFile << "</Placemark>\n";

    if (bIsLinear || bIsHyperbolic || bIsCubic)
    {
        // Get the first and last position waypoints
        const DOMElement *firstPositionWaypointElement = dynamic_cast<const DOMElement *>(positionWaypointList->item(0));
        const DOMElement *lastPositionWaypointElement = dynamic_cast<const DOMElement *>(positionWaypointList->item(positionWaypointList->getLength() - 1));

        // Extract the start and end coordinates
        std::string startCoordinates = getCoordinatesFromPositionWaypoint(firstPositionWaypointElement, referencedLatitude, referencedLongitude, referenceAltitude);
        std::string endCoordinates = getCoordinatesFromPositionWaypoint(lastPositionWaypointElement, referencedLatitude, referencedLongitude, referenceAltitude);

        // Start point placemark
        kmlFile << "<Placemark>\n";
        kmlFile << "    <name>Start: " << XMLString::transcode(element->getAttribute(XMLString::transcode("name"))) << "</name>\n";
        kmlFile << "    <styleUrl>#target</styleUrl>\n"; // Replace with your desired style URL for the start icon
        kmlFile << "    <Point>\n";
        kmlFile << "        <coordinates>" << startCoordinates << "</coordinates>\n";
        if (dAltitudeAboveGround > 0)
        {
            kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
            kmlFile << "        <extrude>1</extrude>\n";
        }
        else
        {
            kmlFile << "        <altitudeMode>clampToGround</altitudeMode>\n";
        }
        kmlFile << "    </Point>\n";
        kmlFile << "</Placemark>\n";

        // End point placemark
        kmlFile << "<Placemark>\n";
        kmlFile << "    <name>End: " << XMLString::transcode(element->getAttribute(XMLString::transcode("name"))) << "</name>\n";
        kmlFile << "    <styleUrl>#target</styleUrl>\n"; // Replace with your desired style URL for the end icon
        kmlFile << "    <Point>\n";
        kmlFile << "        <coordinates>" << endCoordinates << "</coordinates>\n";
        if (dAltitudeAboveGround > 0)
        {
            kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
            kmlFile << "        <extrude>1</extrude>\n";
        }
        else
        {
            kmlFile << "        <altitudeMode>clampToGround</altitudeMode>\n";
        }
        kmlFile << "    </Point>\n";
        kmlFile << "</Placemark>\n";
    }

    return true;
}

// Function to traverse the DOMNode by recursively calling itself and processElement()
void traverseDOMNode(const DOMNode *node, std::ofstream &kmlFile, double referencedLatitude, double referencedLongitude, double referenceAltitude, DOMDocument *document)
{
    if (node->getNodeType() == DOMNode::ELEMENT_NODE)
    {
        const DOMElement *element = dynamic_cast<const DOMElement *>(node);
        if (!processPlatformElement(element, kmlFile, referencedLatitude, referencedLongitude, referenceAltitude, document))
        {
            // Check if the element is a platform
            std::string tagNameStr = xercesc::XMLString::transcode(element->getTagName());
            std::cout << "ERROR: Failed to processing: " + tagNameStr << std::endl;
        }
    }

    for (DOMNode *child = node->getFirstChild(); child != nullptr; child = child->getNextSibling())
    {
        traverseDOMNode(child, kmlFile, referencedLatitude, referencedLongitude, referenceAltitude, document);
    }
}

// Main function
int main(int argc, char *argv[])
{

    // double alpha = 1;
    // double beta = 2;
    // double gamma = 3.6;

    // double angle_3dB_drop_deg = find_3db_drop_angle(alpha, beta, gamma);
    // std::cout << "3dB drop occurs at: " << angle_3dB_drop_deg << " degrees" << std::endl;

    if (argc > 3 && argc < 6)
    {
        std::cerr << "Usage: " << argv[0] << " <input XML file> <output KML file> [<referencedLatitude> <referencedLongitude> <referenceAltitude>]" << std::endl;
        return 1;
    }

    // Setting default geographical and altitude coordinates
    double referencedLatitude = -33.9545;
    double referencedLongitude = 18.4563;
    double referenceAltitude = 0;

    // Update file_path with command line argument
    string file_path = argv[1];
    // Setting mode to evironment variable from command line
    string output_file = argv[2];
    // Setting georgraphical coordinates to command line input
    if (argc == 6)
    {
        try
        {
            referencedLatitude = std::stod(argv[3]);
            referencedLongitude = std::stod(argv[4]);
            referenceAltitude = std::stod(argv[5]);
        }
        catch (const std::invalid_argument &e)
        {
            std::cerr << "Error: Invalid argument. Please provide valid numbers for referencedLatitude, referencedLongitude, and referenceAltitude.\n";
            return 1;
        }
        catch (const std::out_of_range &e)
        {
            std::cerr << "Error: Out of range. Please provide valid numbers for referencedLatitude, referencedLongitude, and referenceAltitude.\n";
            return 1;
        }
    }

    try
    {

        // Initializing Xerces-C++ library:
        XMLPlatformUtils::Initialize();

        // A XercesDOMParser object is set along with its features:
        XercesDOMParser parser;

        // Error handler configuration
        ErrorHandler *errorHandler = (ErrorHandler *)new HandlerBase();
        parser.setErrorHandler(errorHandler);

        // Disables validation during parsing.
        parser.setValidationScheme(XercesDOMParser::Val_Never);

        // Namespace set to false
        parser.setDoNamespaces(false);

        // Validation against schema set to false
        parser.setDoSchema(false);

        parser.setLoadExternalDTD(false);

        // Use file_path from command line argument
        parser.parse(file_path.c_str());

        // Creating DOMDocument and checking if document pointer is valid
        DOMDocument *document = parser.getDocument();
        if (!document)
        {
            std::cerr << "Error: document not found" << std::endl;
            XMLPlatformUtils::Terminate();
            return 1;
        }

        // Creating rootElement and checking if rootElement pointer is valid
        DOMElement *rootElement = document->getDocumentElement();
        if (!rootElement)
        {
            std::cerr << "Error: root element not found" << std::endl;
            XMLPlatformUtils::Terminate();
            return 1;
        }

        std::ofstream kmlFile(output_file.c_str());
        if (!kmlFile.is_open())
        {
            std::cerr << "Error opening output KML file" << std::endl;
            XMLPlatformUtils::Terminate();
            return 1;
        }

        // Write the KML header
        kmlFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        kmlFile << "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\">\n";
        kmlFile << "<Document>\n";
        kmlFile << "<name>" << file_path << "</name>\n";

        // KML styles appended to document
        kmlFile << "<Style id=\"receiver\">\n";
        kmlFile << "  <IconStyle>\n";
        kmlFile << "    <Icon>\n";
        kmlFile << "      <href>https://cdn-icons-png.flaticon.com/512/645/645436.png</href>\n";
        kmlFile << "    </Icon>\n";
        kmlFile << "  </IconStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"transmitter\">\n";
        kmlFile << "  <IconStyle>\n";
        kmlFile << "    <Icon>\n";
        kmlFile << "      <href>https://cdn-icons-png.flaticon.com/128/224/224666.png</href>\n";
        kmlFile << "    </Icon>\n";
        kmlFile << "  </IconStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"target\">\n";
        kmlFile << "  <IconStyle>\n";
        kmlFile << "    <Icon>\n";
        kmlFile << "      <href>https://upload.wikimedia.org/wikipedia/commons/thumb/a/ad/Target_red_dot1.svg/1200px-Target_red_dot1.svg.png</href>\n";
        kmlFile << "    </Icon>\n";
        kmlFile << "  </IconStyle>\n";
        kmlFile << "  <LineStyle>\n";
        kmlFile << "    <width>2</width>\n";
        kmlFile << "  </LineStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"translucentPolygon\">\n";
        kmlFile << "    <LineStyle>\n";
        kmlFile << "        <color>ff0000ff</color>\n";
        kmlFile << "        <width>2</width>\n";
        kmlFile << "    </LineStyle>\n";
        kmlFile << "    <PolyStyle>\n";
        kmlFile << "        <color>00ffffff</color> <!-- RGBA: 50% transparent white --> \n";
        kmlFile << "     </PolyStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"arrowStyle\">\n";
        kmlFile << "    <IconStyle>\n";
        kmlFile << "        <Icon>\n";
        kmlFile << "            <href>http://maps.google.com/mapfiles/kml/shapes/arrow.png</href>\n";
        kmlFile << "        </Icon>\n";
        kmlFile << "        <scale>0.5</scale>\n";
        kmlFile << "    </IconStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"lineStyle\">\n";
        kmlFile << "    <LineStyle>\n";
        kmlFile << "        <color>ff0000ff</color>\n";
        kmlFile << "        <width>2</width>\n";
        kmlFile << "     </LineStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"lineStyleBlue\">\n";
        kmlFile << "    <LineStyle>\n";
        kmlFile << "        <color>ffff0000</color>\n";
        kmlFile << "        <width>2</width>\n";
        kmlFile << "     </LineStyle>\n";
        kmlFile << "</Style>\n";

        // Folder element appended
        kmlFile << "<Folder>\n";
        kmlFile << "  <name>Reference Coordinate</name>\n";
        kmlFile << "  <description>Placemarks for various elements in the FERSXML file. All Placemarks are situated relative to this reference point.</description>\n";

        // Add the LookAt element with given values
        kmlFile << "  <LookAt>\n";
        kmlFile << "    <dLongitude>" << referencedLongitude << "</dLongitude>\n";
        kmlFile << "    <dLatitude>" << referencedLatitude << "</dLatitude>\n";
        kmlFile << "    <altitude>" << referenceAltitude << "</altitude>\n";
        kmlFile << "    <heading>-148.4122922628044</heading>\n";
        kmlFile << "    <tilt>40.5575073395506</tilt>\n";
        kmlFile << "    <range>10000</range>\n";
        kmlFile << "  </LookAt>\n";

        // Traverse DOMNode and output extracted FERSXML data:
        traverseDOMNode(rootElement, kmlFile, referencedLatitude, referencedLongitude, referenceAltitude, document);

        // Close the Folder and Document elements
        kmlFile << "</Folder>\n";
        kmlFile << "</Document>\n";
        kmlFile << "</kml>\n";

        kmlFile.close();

        delete errorHandler; // Clean up the error handler
    }
    catch (const XMLException &e)
    {
        cerr << "Error initializing Xerces-C++: " << XMLString::transcode(e.getMessage()) << endl;
        return 1;
    }
    catch (const DOMException &e)
    {
        cerr << "Error parsing XML: " << XMLString::transcode(e.getMessage()) << endl;
        return 1;
    }
    catch (const SAXException &e)
    {
        cerr << "Error parsing XML: " << XMLString::transcode(e.getMessage()) << endl;
        return 1;
    }
    catch (...)
    {
        cerr << "Unknown error occurred while parsing XML." << endl;
        return 1;
    }

    XMLPlatformUtils::Terminate();
    return 0;
}
