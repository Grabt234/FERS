import xml.etree.ElementTree as ET

class SignalGenerator:

    def generate_chirp(ramp_period, sample_rate, bandwidth, init_freq=0, tau=0, phi=0):
        """Generate an analytic baseband chirp.

        Args:
            ramp_period (float): Ramp period of chirp (seconds)
            sample_rate (uint): Sample rate used to generate chirp (Hz)
            bandwidth (uint): Bandwidth of chirp (Hz)
            init_freq (int, optional): Start frequency (Hz). Defaults to 0.
            tau (int, optional): Sweep offset (s). Defaults to 0.
            phi (int, optional): Phase offset (rad). Defaults to 0.

        Returns:
            np.array: chirp data array
        """

        import numpy as np

        ns_chirp = int(np.ceil(ramp_period*sample_rate))
        t_chirp = np.linspace(-ramp_period/2, ramp_period/2, ns_chirp, endpoint=False)


        return np.exp(1.j*(np.pi*bandwidth/(2*max(t_chirp))*pow((t_chirp - tau), 2) + 2*np.pi*init_freq*(t_chirp - tau) + phi))
    
    def write_hdf5(IQ_data, file_name):
        """Write IQ data to an HDF5 file.

        Args:
            IQ_data (np.array): numpy array of complex data
            filename (string): name file written file
        """

        import os as os

        _, file_extension = os.path.splitext(file_name)
        if file_extension != ".h5":
            file_name = file_name + ".h5"

        import h5py
        import numpy as np

        h5 = h5py.File(file_name, 'w')
        h5.create_dataset('/I/value', data=np.real(IQ_data))
        h5.create_dataset('/Q/value', data=np.imag(IQ_data))
        h5.close()


    def read_hdf5(filename):
        """ Read IQ data from an HDF5 file.

        Args:
            filename (string): path to .h5 file to read in 

        Returns:
            np.array: Array of complex IQ data
        """

        import os
        import h5py
        import numpy as np

        if (os.path.exists(filename) == False):
            print("HDF5 file not found. Please check the path.")
            exit()
        
        h5 = h5py.File(filename, 'r')

        dataset_list = list(h5.keys())

        # read attributes
        # attribute_list = h5[dataset_list[0]].attrs.keys()
        # for attr in attribute_list:
            # print(attr, h5[dataset_list[0]].attrs[attr])

        scale = np.float64(h5[dataset_list[0]].attrs['fullscale'])
        # rate = np.float64(h5[dataset_list[0]].attrs['rate'])
        # time = np.float64(h5[dataset_list[0]].attrs['time'])

        n_pulses = int(np.floor(np.size(dataset_list)/2))
        ns_pulse = int(np.size(h5[dataset_list[0]]))

        i_matrix = np.zeros((n_pulses, ns_pulse), dtype='float64')
        q_matrix = np.zeros((n_pulses, ns_pulse), dtype='float64')

        for i in range(0, n_pulses):
            i_matrix[i, :] = np.array(h5[dataset_list[2*i + 0]], dtype='float64')
            q_matrix[i, :] = np.array(h5[dataset_list[2*i + 1]], dtype='float64')

        dataset = np.array(i_matrix + 1j*q_matrix).astype('complex128')

        dataset *= scale
        
        return dataset
    
    def read_wav(file_name):
        """
        Reads a WAV audio file and returns the audio data and sample rate.

        Args:
            file_name (str): The path to the WAV file.

        Returns:
            tuple: A tuple containing the audio data (numpy array) and the sample rate (int).
        """
        
        import wave
        
        with wave.open(file_name, 'r') as wav:
            # Get the number of channels, sample width, frame rate, and number of frames
            num_channels = wav.getnchannels()
            sample_width = wav.getsampwidth()
            sample_rate = wav.getframerate()
            num_frames = wav.getnframes()

            # Read the audio data
            audio_data = wav.readframes(num_frames)

        # Convert the audio data to a numpy array
        if sample_width == 1:
            audio_data = numpy.frombuffer(audio_data, dtype=numpy.int8)
        elif sample_width == 2:
            audio_data = numpy.frombuffer(audio_data, dtype=numpy.int16)
        elif sample_width == 4:
            audio_data = numpy.frombuffer(audio_data, dtype=numpy.int32)
        else:
            raise ValueError(f"Unsupported sample width: {sample_width}")

        # Reshape the audio data to have the correct number of channels
        audio_data = audio_data.reshape(-1, num_channels)

        return audio_data, sample_rate
        
class SimulationConfiguration:
    def __init__(self, name):
        """Object to store simulation description

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
        self.dtd = ""

    def set_simulaion_parameters(self, start_time, end_time, propogation_speed, interprate, sample_rate, adc_bits):
        """ Set the simulation parameters of the scenario

        Args:
            start_time (string): time when recording of simulation starts
            end_time (string): End time of simulation
            propogation_speed (string): Speed of propogation of waves
            interprate (string):  
            sample_rate (string): 
            adc_seed (string): ADC bits used to simulate recorded signals
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

        param = ET.SubElement(self.parameters, 'adc_bits')
        param.text = str(adc_bits)



    def define_pulse(self, name, filename, length, power, carrier):
        """Define a pulse

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
        timing = ET.SubElement(self.simulation, 'timing', name=name, synconpulse=sync_on_pulse)
        
        ET.SubElement(timing, 'frequency').text = str(frequency)
       
        ET.SubElement(timing, 'freq_offset').text = str(freq_offset)
        ET.SubElement(timing, 'phase_offset').text = str(phase_offset)

        if (freq_offset != "undefined"):
            ET.SubElement(timing, 'random_freq_offset').text = str(std_dev_freq_offset)
        else:
            ET.SubElement(timing, 'random_freq_offset').text = str(0)
        
        if (phase_offset != "undefined"):
            ET.SubElement(timing, 'random_phase_offset').text = str(std_dev_phase_offset)
        else:
            ET.SubElement(timing, 'random_phase_offset').text = str(0)
        
        if (alpha != "undefined" and weight != "undefined"):
            noise = ET.SubElement(timing, 'noise_entry')
            ET.SubElement(noise, 'alpha1').text = str(alpha)
            ET.SubElement(noise, 'weight').text =  str(weight)

        self.pulses.append(timing)


    def define_isotropic_antenna(self, name, efficiency):
        """Define antenna

        Args:
            name (string): name of the antenna
            efficiency (string): efficiency of antenna (0-1)
        """
        
        ant = ET.SubElement(self.simulation, "antenna", pattern='isotropic', name=name)
        ET.SubElement(ant, 'efficiency').text = str(efficiency)


    def define_parabolic_antenna(self, name, efficiency, diameter):
        """Define antenna

        Args:
            name (string): name of the antenna
            efficiency (string): efficiency of antenna
            diameter (string): diameter of antenna ine meters
        """

        ant = ET.SubElement(self.simulation, "antenna", pattern='parabolic', name=name)
        ET.SubElement(ant, 'diameter').text = str(diameter)
        ET.SubElement(ant, 'efficiency').text = str(efficiency)

    def define_sync_antenna(self, name, efficiency, alpha, beta, gamma):
        """Define antenna

        Args:
            TODO: dont know how this works :)
        """

        ant = ET.SubElement(self.simulation, "antenna", pattern='sinc', name=name)
        ET.SubElement(ant, 'alpha').text = str(alpha)
        ET.SubElement(ant, 'beta').text = str(beta)
        ET.SubElement(ant, 'gamma').text = str(gamma)
        ET.SubElement(ant, 'efficiency').text = str(efficiency)
    
    def create_static_isotropic_target_platform(self, platform_name, target_name, rcs, x, y, z):
        """Add a static isotropic target to simulation

        Args:
            platform_name (string): platform name
            target_name (string): target name
            rcs (string): RCS of target
            x (string): x position of platform
            y (string): y position of platform
            z (string): Z position of platform
        """

        platform = ET.SubElement(self.simulation, "platform", name=platform_name)
        
        motion_path = ET.SubElement(platform, "motionpath",interpolation="linear")
        position_waypoint = ET.SubElement(motion_path, "positionwaypoint")
        ET.SubElement(position_waypoint, 'x').text = str(x)
        ET.SubElement(position_waypoint, 'y').text = str(y)
        ET.SubElement(position_waypoint, 'altitude').text = str(z)
        ET.SubElement(position_waypoint, 'time').text = str(0.000001)

        fixedrotation = ET.SubElement(platform, "fixedrotation")
        ET.SubElement(fixedrotation, 'startazimuth').text = str(0)
        ET.SubElement(fixedrotation, 'startelevation').text = str(0)
        ET.SubElement(fixedrotation, 'azimuthrate').text = str(0)
        ET.SubElement(fixedrotation, 'elevationrate').text = str(0)

        target = ET.SubElement(platform, "target", name=target_name)
        rcs_element = ET.SubElement(target, "rcs", type="isotropic")
        ET.SubElement(rcs_element, "value").text = str(rcs)

    def __add_waypoint(self, path, x, y, z, t):
        """Add a waypoint in the fersxml

        Args:
            path (string): path to xml document to write
            x (int): X position at a point in time
            y (int): Y position at a point in time
            z (int): Z position at a point in time
            t (int): Time point where these positions occur
        """

        point = ET.SubElement(path, "positionwaypoint")

        t_x = ET.SubElement(point, 'x')
        t_x.text = str(x)

        t_y = ET.SubElement(point, 'y')
        t_y.text = str(y)

        t_a = ET.SubElement(point, 'altitude')
        t_a.text = str(z)

        t_t = ET.SubElement(point, 'time')
        t_t.text = str(t)

    def create_moving_isotropic_target_platform(self, platform_name, target_name, rcs, x, y, z, t):
        """Define an isotropic radiator

        Args:
            name (string): name of target
            rcs (string): rcs of target
        """

        platform = ET.SubElement(self.simulation, "platform", name=platform_name)
        
        motion_path = ET.SubElement(platform, "motionpath",interpolation="linear")
        for i in range(len(x)):
            self.__add_waypoint(motion_path, x[i], y[i], z[i], t[i])

        fixedrotation = ET.SubElement(platform, "fixedrotation")
        ET.SubElement(fixedrotation, 'startazimuth').text = str(0)
        ET.SubElement(fixedrotation, 'startelevation').text = str(0)
        ET.SubElement(fixedrotation, 'azimuthrate').text = str(0)
        ET.SubElement(fixedrotation, 'elevationrate').text = str(0)

        target = ET.SubElement(platform, "target", name=target_name)
        rcs_element = ET.SubElement(target, "rcs", type="isotropic")
        ET.SubElement(rcs_element, "value").text = str(rcs)

    def create_static_receiver_platform(self, platform_name, receiver_name, antenna, timing, nodirect, nopropagationloss, window_skip, window_length, prf, noise_temp, x, y, z, boresight_azimuth, boresight_elevation):
        """Create a static reciever platform

        Args:
            platform_name (string): platform name used by fers
            receiver_name (string): receiver name used by fers
            # antenna (string): name of timing source as previously defined
            # timing (string): name of timing source as previously defined
            # nodirect (string): "true" or "false" if direct path signal is ignored
            # nopropagationloss (string): "true" or "false" causes propagation loss to be ignored
            # window_skip (string): Time to skip after start of pulse before starting receiving (seconds)
            # window_length (string): Length of the range gate (seconds)
            # prf (string): 1/window_length
            # noise_temp (string): Noise temperature of the receiver
            x (string): X Position of transmitter (meter)
            y (string): Y Position of transmitter (meter)
            z (string): Z Position of transmitter (meter)
            t (string): T in instance (seconds) at each previously defined position
            boresight_azimuth (string): Boresight azimuth of antenna (can be 0 if isotropic)
            boresight_elevation (string): Boresight elevation of antenna (can be 0 if isotropic)
        """

        platform = ET.SubElement(self.simulation, "platform", name=platform_name)

        motion_path = ET.SubElement(platform, "motionpath",interpolation="linear")
        position_waypoint = ET.SubElement(motion_path, "positionwaypoint")
        ET.SubElement(position_waypoint, 'x').text = str(x)
        ET.SubElement(position_waypoint, 'y').text = str(y)
        ET.SubElement(position_waypoint, 'altitude').text = str(z)
        ET.SubElement(position_waypoint, 'time').text = str(0.000001)

        fixedrotation = ET.SubElement(platform, "fixedrotation")
        ET.SubElement(fixedrotation, 'startazimuth').text = str(boresight_azimuth)
        ET.SubElement(fixedrotation, 'startelevation').text = str(boresight_elevation)
        ET.SubElement(fixedrotation, 'azimuthrate').text = str(0)
        ET.SubElement(fixedrotation, 'elevationrate').text = str(0)

        receiver = ET.SubElement(platform, "receiver", name=receiver_name, antenna=antenna, timing=timing, nodirect=nodirect, nopropagationloss=nopropagationloss)
        ET.SubElement(receiver, 'window_skip').text = str(window_skip)
        ET.SubElement(receiver, 'window_length').text = str(window_length)
        ET.SubElement(receiver, 'prf').text = str(prf)
        ET.SubElement(receiver, 'noise_temp').text = str(noise_temp)

    def create_moving_receiver_platform(self, platform_name, receiver_name, antenna, timing, nodirect, nopropagationloss, window_skip, window_length, prf, noise_temp, x, y, z,t, boresight_azimuth, boresight_elevation):
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

        platform = ET.SubElement(self.simulation, "platform", name=platform_name)

        motion_path = ET.SubElement(platform, "motionpath",interpolation="linear")
        for i in range(len(x)):
            self.__add_waypoint(motion_path, x[i], y[i], z[i], t[i])

        fixedrotation = ET.SubElement(platform, "fixedrotation")
        ET.SubElement(fixedrotation, 'startazimuth').text = str(boresight_azimuth)
        ET.SubElement(fixedrotation, 'startelevation').text = str(boresight_elevation)
        ET.SubElement(fixedrotation, 'azimuthrate').text = str(0)
        ET.SubElement(fixedrotation, 'elevationrate').text = str(0)

        receiver = ET.SubElement(platform, "receiver", name=receiver_name, antenna=antenna, timing=timing, nodirect=nodirect, nopropagationloss=nopropagationloss)
        ET.SubElement(receiver, 'window_skip').text = str(window_skip)
        ET.SubElement(receiver, 'window_length').text = str(window_length)
        ET.SubElement(receiver, 'prf').text = str(prf)
        ET.SubElement(receiver, 'noise_temp').text = str(noise_temp)

    def create_static_transmitter_platform(self, platform_name, transmitter_name,pulse_name, antenna, timing, type, prf, x, y, z, boresight_azimuth, boresight_elevation):
        """Create a static transmitter platform

        Args:
            platform_name (string): platform name used by fers
            transmitter_name (string): transmitter name used by fers
            pulse_name (string): Pulse name previously defined in simulation setup and used by transmitter 
            antenna (string): Antenna name previously defined in simulation setup and used by transmitter
            timing (string): Timing source name previously defined in simulation setup and used by transmitter
            type (string): "continous" or "pulsed" transmission type
            prf (string): If pulsed define the PRF of the transmission 
            x (string): X Position of transmitter (meter)
            y (string): Y Position of transmitter (meter)
            z (string): Z Position of transmitter (meter)
            t (string): T in instance (seconds) at each previously defined position
            boresight_azimuth (string): Boresight azimuth of antenna (can be 0 if isotropic)
            boresight_elevation (string): Boresight elevation of antenna (can be 0 if isotropic)
        """

        platform = ET.SubElement(self.simulation, "platform", name=platform_name)

        motion_path = ET.SubElement(platform, "motionpath",interpolation="linear")
        position_waypoint = ET.SubElement(motion_path, "positionwaypoint")
        ET.SubElement(position_waypoint, 'x').text = str(x)
        ET.SubElement(position_waypoint, 'y').text = str(y)
        ET.SubElement(position_waypoint, 'altitude').text = str(z)
        ET.SubElement(position_waypoint, 'time').text = str(0.000001)

        fixedrotation = ET.SubElement(platform, "fixedrotation")
        ET.SubElement(fixedrotation, 'startazimuth').text = str(boresight_azimuth)
        ET.SubElement(fixedrotation, 'startelevation').text = str(boresight_elevation)
        ET.SubElement(fixedrotation, 'azimuthrate').text = str(0)
        ET.SubElement(fixedrotation, 'elevationrate').text = str(0)

        receiver = ET.SubElement(platform, "transmitter", name=transmitter_name, type=type, antenna=antenna, timing=timing, pulse=pulse_name)
        ET.SubElement(receiver, 'prf').text = str(prf)

    def create_moving_transmitter_platform(self, platform_name, transmitter_name,pulse_name, antenna, timing, type, prf, x, y, z,t, boresight_azimuth, boresight_elevation):
        """Create a moving transmitter platform

        Args:
            platform_name (string): platform name used by fers
            transmitter_name (string): transmitter name used by fers
            pulse_name (string): Pulse name previously defined in simulation setup and used by transmitter 
            antenna (string): Antenna name previously defined in simulation setup and used by transmitter
            timing (string): Timing source name previously defined in simulation setup and used by transmitter
            type (string): "continous" or "pulsed" transmission type
            prf (string): If pulsed define the PRF of the transmission 
            x (list[int]): X Positions of transmitter (meter)
            y (list[int]): Y Positions of transmitter (meter)
            z (list[int]): Z Positions of transmitter (meter)
            t (list[int]): T in instances (seconds) at each previously defined position
            boresight_azimuth (string): Boresight azimuth of antenna (can be 0 if isotropic)
            boresight_elevation (string): Boresight elevation of antenna (can be 0 if isotropic)
        """

        platform = ET.SubElement(self.simulation, "platform", name=platform_name)

        motion_path = ET.SubElement(platform, "motionpath",interpolation="linear")
        for i in range(len(x)):
            self.__add_waypoint(motion_path, x[i], y[i], z[i], t[i])

        fixedrotation = ET.SubElement(platform, "fixedrotation")
        ET.SubElement(fixedrotation, 'startazimuth').text = str(boresight_azimuth)
        ET.SubElement(fixedrotation, 'startelevation').text = str(boresight_elevation)
        ET.SubElement(fixedrotation, 'azimuthrate').text = str(0)
        ET.SubElement(fixedrotation, 'elevationrate').text = str(0)

        receiver = ET.SubElement(platform, "transmitter", name=transmitter_name, type=type, antenna=antenna, timing=timing, pulse=pulse_name)
        ET.SubElement(receiver, 'prf').text = str(prf)

    def set_export_options(self, xml="true", csv="true", binary="false", csvbinary="false"):
        """Set what shall be exported from simulation

        Args:
            xml (str, optional): export xml. Defaults to "true".
            csv (str, optional): export csv. Defaults to "true".
            binary (str, optional): export binary. Defaults to "false".
            csvbinary (str, optional): export csv binary. Defaults to "false".
        """
        export = ET.SubElement(self.parameters, 'export', xml=xml, csv=csv, binary=binary, csvbinary=csvbinary)

    def to_xml_string(self):
        """get the simulation as an xml string

        Returns:
            string: returns simulation as an xml string
        """
        return ET.tostring(self.simulation, encoding='utf-8')
    
    def define_dtd_file(self, name):
        """Define DTD for xml using inpit

        Args:
            name (string): DTD definition file for xml
        """
        self.dtd = doctype_declaration = "<!DOCTYPE simulation SYSTEM \""+name+"\" >"

    def write_to_file(self, file_path):
        """Write the simulation to a fersxml file

        Args:
            file_path (string): path to which to write fersxml
        """
        
        # Parse the XML string
        import xml.dom.minidom
        import xml.etree.ElementTree as ET

        xml_string = ET.tostring(self.simulation, encoding='utf-8')
        
        dom = xml.dom.minidom.parseString(xml_string)
        pretty_xml = dom.toprettyxml(indent="  ")

        # Split the string into lines
        lines = pretty_xml.split('\n')
        # Remove the first line
        lines.pop(0)
        # Join the remaining lines back together
        pretty_xml = '\n'.join(lines)

        with open(file_path, "wb") as f:
            f.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n".encode('utf-8'))
            f.write(self.dtd.encode('utf-8'))
            f.write(pretty_xml.encode('utf-8'))

# only contains single reflector
class SimpleStaticIsotropicTargetPlatform:
    def __init__(self, platform_name, target_name, x, y, z, rcs):
        self.platform_name = platform_name
        self.target_name = target_name
        self.x = x
        self.y = y
        self.z = z
        self.rcs = rcs

class SimpleMovingIsotropicTargetPlatform:
    def __init__(self, platform_name, target_name, x, y, z, t, rcs):
        self.platform_name = platform_name
        self.target_name = target_name
        self.x = x
        self.y = y
        self.z = z
        self.t = t
        self.rcs = rcs