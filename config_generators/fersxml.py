import xml.etree.ElementTree as ET

class SimulationControl:
    def __init__(self, name):
        """_summary_

        Args:
            name (string): name of the fers simulation
        """
        self.simulation = ET.Element('simulation', name=name)
        self.parameters = ET.SubElement(self.simulation, 'parameters')
        self.pulses = []
        self.waveforms = []
        self.timings = []
        self.antennas = []
        self.platforms = []
        self.includes = []

    def set_simulaion_parameters(self, start_time, end_time, propogation_speed, interprate, sample_rate, random_seed, adc_bits):
        """ Set the simulation parameters of the scenario

        Args:
            start_time (string): time when recording of simulation starts
            end_time (string): End time of simulation
            propogation_speed (string): Start time of simulation
            interprate (string):  Position interpolation rate for CW
            sample_rate (string): Override the rendering sample rate with the specified one (Hz)
            random_seed (string): Random seed for noise and jitter (positive integer). If this is not specified, time() is used
            adc_seed (string): Formats to export
        """

        param = ET.SubElement(self.parameters, 'starttime')
        param.text = str(start_time)

        param = ET.SubElement(self.parameters, 'endtime')
        param.text = str(end_time)

        param = ET.SubElement(self.parameters, 'c')
        param.text = str(propogation_speed)

        param = ET.SubElement(self.parameters, 'interprate')
        param.text = str(interprate)

        param = ET.SubElement(self.parameters, 'rate')
        param.text = str(sample_rate)

        param = ET.SubElement(self.parameters, 'randomseed')
        param.text = str(random_seed)

        param = ET.SubElement(self.parameters, 'adcbits')
        param.text = str(adc_bits)



    def define_pulse(self, name, filename, length, power, carrier):
        """Define a pulse

            note: rate parameter is left default as its function is unclear

        Args:
            name (string): Name of the pulse to be references in other parts of simulation
            filename (string): File name to be transmitted (HDF5)
            length (string): Length of pulse in seconds (I think dead time + pulse period (i.e 1/hop rate))
            power (string): Power in Watts
            carrier (string): Center frequency in Hz
        """
        pulse = ET.SubElement(self.simulation, 'pulse', name=name, type="file", filename=filename)
        
        ET.SubElement(pulse, 'length').text = str(length)
        ET.SubElement(pulse, 'power').text = str(power)
        ET.SubElement(pulse, 'carrier').text = str(carrier)
        
        self.pulses.append(pulse)

    # Add similar methods for waveform, timing, antenna, platform, and include elements

    def set_export_options(self, xml="true", csv="true", binary="false", csvbinary="false"):
        export = ET.SubElement(self.parameters, 'export', xml=xml, csv=csv, binary=binary, csvbinary=csvbinary)

    def to_xml_string(self):
        return ET.tostring(self.simulation, encoding='unicode')

    def write_to_file(self, file_path):
        tree = ET.ElementTree(self.simulation)
        tree.write(file_path, encoding='unicode', xml_declaration=True)