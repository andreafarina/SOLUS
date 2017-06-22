function data = f_read_sdt_01(fn) 
% hw 28.01.12
% read sdt files - for fgSPAD scanner campaign only 
% (2 data blocks)
% file header has 1365 bytes
% 22-byte data block headers before every data block
% reads data only, does not yet handle all other info in the files
% so far for these specific 2-block files only
% structure derived from SPC_data_file_structure.h 23/06/11
% from TCSPC package 3.5.1
% 01.07.13 display suppressed

fid = fopen(fn,'r');

% hope with this replacement it will also work on 64 bit machine
synon_ulong = 'ubit32'; %'ulong';
synon_long = 'bit32'; %'long';
  
%% 1 - File header

% typedef struct  {
%    short    revision;  // software revision & module identification
%                        //   lowest bits 0-3 -   software revision ( >= 12(decimal))
%                        //   bits 11-4   - BH module type
%                        //      meaning of this field values (hex):
%                        //        0x20 -SPC-130, 0x21 - SPC-600, 0x22 - SPC-630,
%                        //        0x23 -SPC-700, 0x24 - SPC-730, 0x25 - SPC-830,
%                        //        0x26 -SPC-140, 0x27 - SPC-930, 0x28 - SPC-150,
%                        //        0x29 -DPC-230, 0x2a - SPC-130EM
%                        //   highest bits 15-12 - module subtype - not used yet
%    long     info_offs; // offset of the info part which contains general text 
%                        //   information (Title, date, time, contents etc.)
%    short    info_length;  // length of the info part
%    long     setup_offs;   // offset of the setup text data 
%                // (system parameters, display parameters, trace parameters etc.)
%    short    setup_length;  // length of the setup data
%    long     data_block_offs;   // offset of the first data block 
%    short    no_of_data_blocks; // no_of_data_blocks valid only when in 0 .. 0x7ffe range,
%                                // if equal to 0x7fff  the  field 'reserved1' contains 
%                                //     valid no_of_data_blocks
%    long     data_block_length;     // length of the longest block in the file
%    long     meas_desc_block_offs;  // offset to 1st. measurement description block 
%                                    //   (system parameters connected to data blocks)
%    short    no_of_meas_desc_blocks;  // number of measurement description blocks
%    short    meas_desc_block_length;  // length of the measurement description blocks
%    unsigned short    header_valid;   // valid: 0x5555, not valid: 0x1111
%    unsigned long     reserved1;      // reserved1 now contains no_of_data_blocks
%                                      // reserved2 now contains length (in int words) of data block extension
%    unsigned short    reserved2;
%    unsigned short    chksum;            // checksum of file header
%    }bhfile_header;

  revision      = fread(fid, 1, 'short');
  info_offs     = fread(fid, 1, synon_long);
  info_length   = fread(fid, 1, 'short');
  setup_offs    = fread(fid, 1, synon_long);
  setup_length  = fread(fid, 1, 'short');
  data_block_offs = fread(fid, 1, synon_long);
  no_of_data_blocks = fread(fid, 1, 'short');
  data_block_length  = fread(fid, 1, synon_long);
  meas_desc_block_offs  = fread(fid, 1, synon_long);
  no_of_meas_desc_blocks = fread(fid, 1, 'short');
  meas_desc_block_length = fread(fid, 1, 'short');
  header_valid = fread(fid, 1, 'ushort');
  no_of_data_blocks_reserved1 = fread(fid, 1, synon_ulong);
  reserved2 = fread(fid, 1, 'ushort'); % obviously not
  chksum = fread(fid, 1, 'ushort');

%% 2 - File Info

infotext = char(fread(fid, info_length, 'char')');

%% 3 - Setup
% seems to be not in use in the files measured

% The setup block contains all the system parameters, display parameters, trace parameters etc 
%   stored in ASCII.
% From the software version 8.1  setup block contains also binary part.
%   
% It is used to set the SPC system (hardware and software) into the same state as 
%   it was in the moment when the data file was stored. 
% The values are stored together with an identifier of the particular parameter.
% This method allows to maintain compatibility between different SPC versions. 
% If a parameter is missing in the setup part, i.e. 
%   if a file from an older software version is loaded, a default value is used 
%    when the file is loaded. 
% For Multi SPC Systems the system parameters section contains subsections 
%    for module parameters which are separate for the individual modules.
% 
% ASCII part of the setup block starts with the line
% *SETUP
% and ends with the line
% *END
% 
%  DI_ parameters in ASCII part of the setup are now obsolete and are present
%  for compatibility reason only - parameters from binary setup are used
% 
% After ASCII part follows the binary part of the setup  ( for the files created 
%     with the software version 8.1 and next ) 
% 
% It starts with the text BIN_PARA_BEGIN:
% Directly after BIN_PARA_BEGIN: text comes:
%   -  1 byte = '0' which is EOS character ( End Of String )
%   -  4-byte value which is equal to the length of binary setup 
%           ( all suqsequent headers and structures up to *END text )
%   -  BHBinHdr bh_bin_hdr       
%   -  SPCBinHdr spc_bin_hdr       
%   -  structures defined in spc_bin_hdr header
% 
% binary setup ends with the line
% *END
% 
% binary setup contains mainly display parameters for different graph windows
%  ( Main graph (2D, 3D), FCS and other (in future) ), 
%  which are used by the main software to restore display parameters
%    of the specific windows
%  It contains also parameters specific for external devices controlled by SPCM software
%    like GVD-120 module, SyncDelay Box or parameters for Zscan with Zeiss AxioObserver
% 
% In order to reach binary part of the setup user should read the whole setup part 
%   to the buffer and then read strings from the buffer 
%     until the string BIN_PARA_BEGIN: is found.

%% 4 - Measurement Description Blocks
% seems to be not in use in the files measured
  
% Each data block can (but need not) have its own system (hardware) parameter set 
%  which can differ from the setup parameters. 

%%  5- Data Blocks

% With the software version 7.0  the data block header was changed to make possible 
%   a higher number of data blocks and a variable block size.
% Each data block can now contain a 'Data Set' i.e. the data of several curves 
%   which were obtained in one measurement. 
% % % The number and the location of  the data blocks is contained in the file header 
% % %   at the beginning of the data file. 
% % % The length of the block is contained in the block header.
% The data block header contains the data block number, 
%   the offset of the data block from the beginning of the file, 
%   the offset to the next data block, 
%   an information about the data in the block (measured block, block loaded from file, etc.), 
%   and a reference to the corresponding measurement description block
% 
% The data block header structure is shown below.
% */
% 
% 
% typedef struct {
%    short          block_no;   // number of the block in the file
%                          // valid only  when in 0 .. 0x7ffe range, otherwise use lblock_no field
%                          // obsolete now, lblock_no contains full block no information                         
%    long           data_offs;       // offset of the data block from the beginning of the file
%    long           next_block_offs; // offset to the data block header of the next data block
%    unsigned short block_type;      // see block_type defines below
%    short          meas_desc_block_no; // Number of the measurement description block 
%                                       //    corresponding to this data block
%    unsigned long  lblock_no;       // long block_no - see remarks below 
%    unsigned long  block_length;    // reserved2 now contains block( set ) length
%    }BHFileBlockHeader;

status = fseek(fid, data_block_offs, 'bof'); % bof = beginning of file

  % typedef struct {
  %    short          block_no;   // number of the block in the file
  %                          // valid only  when in 0 .. 0x7ffe range, otherwise use lblock_no field
  %                          // obsolete now, lblock_no contains full block no information                         
  %    long           data_offs;       // offset of the data block from the beginning of the file
  %    long           next_block_offs; // offset to the data block header of the next data block
  %    unsigned short block_type;      // see block_type defines below

  % 2 - continuos flow measurement ( BIFL )  

  %    short          meas_desc_block_no; // Number of the measurement description block 
  %                                       //    corresponding to this data block
  %    unsigned long  lblock_no;       // long block_no - see remarks below 
  %    unsigned long  block_length;    // reserved2 now contains block( set ) length
  %    }BHFileBlockHeader;

  
seed = uint16(0);
data = repmat(seed, data_block_length / 2, no_of_data_blocks);


%% reading data
% note that this is just for a 2-block file, not yet general

for i = 1:no_of_data_blocks
  block_no = fread(fid, 1, 'short');
  data_offs = fread(fid, 1, synon_long);
  next_block_offs = fread(fid, 1, synon_long);
  block_type = fread(fid, 1, 'ushort');
  meas_desc_block_no = fread(fid, 1, 'short');
  lblock_no = fread(fid, 1, synon_ulong);
  block_length = fread(fid, 1, synon_ulong);
    
  if i == 1
    status = fseek(fid, data_offs, 'bof'); % bof = beginning of file
    hh = fread(fid, data_block_length / 2, 'ushort'); % if length is in bytes
    status = fseek(fid, next_block_offs, 'bof'); % bof = beginning of file
  
  elseif i == 2
    hh = fread(fid, data_block_length / 2, 'ushort'); % if length is in bytes
  end
  data(:,i) = hh;
  
end;

fclose(fid);


