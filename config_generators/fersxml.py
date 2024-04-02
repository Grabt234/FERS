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


    # Add similar methods for waveform, timing, antenna, weight, and include elements
    def define_timing_source(self, name, sync_on_pulse, alpha, weight,frequency, freq_offset, std_dev_freq_offset, phase_offset, std_dev_phase_offset):
        """Define Timing source

        Args:
            name (string): name of the clock source
            sync_on_pulse (string): set to "true" or "false" ( Reset timing error at the start of each pulse)
            alpha (string): Alpha parameters for 1/f^alpha clock model, set to "undefined" if unused
            weight (string)://!< Weights for 1/f^alpha clock model, set to "undefined" if unused
            frequency (string): The nominal oscillator frequency
            freq_offset (string): frequency offset of clock, , set to "undefined" if unused
            std_dev_freq_offset (string): stdev of phase offset, if freq_offset defined, this shall not be written
            phase_offet (string): Phase offset of clock
            std_dev_phase_offset (string): stdev of phase offset, if phase _offset defined, this shall not be written
            
        """
        timing = ET.SubElement(self.simulation, 'timing', name=name, sync_on_pulse=sync_on_pulse)
        
        ET.SubElement(timing, 'frequency').text = str(frequency)

        if (alpha != "undefined" and weight != "undefined"):
            ET.SubElement(timing, 'noise_entry',alpha=alpha, weight=weight)

        
        ET.SubElement(timing, 'freq_offset').text = str(freq_offset)

        if (freq_offset != "undefined"):
            ET.SubElement(timing, 'random_freq_offset').text = str(std_dev_freq_offset)

        ET.SubElement(timing, 'phase_offset').text = str(phase_offset)
        
        if (phase_offset != "undefined"):
            ET.SubElement(timing, 'random_phase_offset').text = str(std_dev_phase_offset)
        
        self.pulses.append(timing)


    def define_isotropic_antenna(self, name, efficiency):
        """Define antenna

        Args:
            name (string): name of the antenna
            efficiency (string): efficiency of antenna (0-1)
        """
        
        antenna = ET.SubElement(self.simulation, 'isotropic', name=name)
        ET.SubElement(antenna, 'efficiency').text = str(efficiency)


    def define_parabolic_antenna(self, name, efficiency, diameter):
        """Define antenna

        Args:
            name (string): name of the antenna
            efficiency (string): efficiency of antenna
            diameter (string): diameter of antenna ine meters
        """

        antenna = ET.SubElement(self.simulation, 'parabolic', name=name)
        ET.SubElement(antenna, 'efficiency').text = str(efficiency)
        ET.SubElement(antenna, 'diameter').text = str(diameter)
    
    def define_isotropic_target(self, name, rcs):
        """Define an isotropic radiator

        Args:
            name (string): name of target
            rcs (string): rcs of target
        """

        target = ET.SubElement(self.simulation, "target", name=name)
        ET.SubElement(target, "rcs", type="isotropic")
        ET.SubElement(target, "value").text = str(rcs)

    def define_receiver(self, name, antenna, timing, nodirect, nopropagationloss, window_skip, window_length, prf, noise_temp):
        """Define a reciever with specified clock and antenna

        Args:
            name (string): name of the receiver
            antenna (string): name of timing source as previously defined
            timing (string): name of timing source as previously defined
            nodirect (string): "true" or "false" if direct path signal is ignored
            nopropagationloss (string): "true" or "false" causes propagation loss to be ignored
            window_skip (string): Time to skip after start of pulse before starting receiving (seconds)
            window_length (string): Length of the range gate (seconds)
            prf (string): 1/window_length
            noise_temp (string): Noise temperature of the receiver
        """

        receiver = ET.SubElement(self.simulation, "receiver", name=name, antenna=antenna, timing=timing, nodirect=nodirect, nopropagationloss=nopropagationloss)
        ET.SubElement(receiver, 'window_skip').text = str(window_skip)
        ET.SubElement(receiver, 'window_length').text = str(window_length)
        ET.SubElement(receiver, 'prf').text = str(prf)
        ET.SubElement(receiver, 'noise_temp').text = str(noise_temp)


    def set_export_options(self, xml="true", csv="true", binary="false", csvbinary="false"):
        export = ET.SubElement(self.parameters, 'export', xml=xml, csv=csv, binary=binary, csvbinary=csvbinary)

    def to_xml_string(self):
        return ET.tostring(self.simulation, encoding='unicode')

    def write_to_file(self, file_path):
        tree = ET.ElementTree(self.simulation)
        tree.write(file_path, xml_declaration=True)