/* pcidev2.h */
# ifndef _PCIDEV_
# define _PCIDEV_

# define bit(n) (1 << (n))
# define bits(pat,n) ((pat) << (n))
# define bitinv(n) (~(1 << (n)))
# define bitsinv(pat,n) (~((pat) << (n)))

enum { MaxPorts=4, PortSize=2, MemWidth=4 };
enum { OFF, ON };
enum { SetIOMode=100, SetModes, SetPortSize, SetFreq, SetTrigger,
       SetOnePort, SetExtraBufs,
       GetPCIConfig=200, GetPlxConfig, GetXilinxConfig,
       GetPortInfo, GetPortAddr, GetBlockInfo,
       ReadMemData, WriteMemData
     };
enum { ConfigWords=64, ConfigSize=ConfigWords*MemWidth,
       PCIConfWords=16, PCIConfSize=PCIConfWords*MemWidth,
       ClassDefs=4, BaseAddrs=6,
       AddrSize=sizeof(void*)
     };
enum XilinxClkPar
     { InternalVarClk=0, ExternalClk=1, Internal10=2, Internal20=3,
       MaxClkSet=10000000 /* Same as Internal10 */
     };
enum XilinxTrigAddr
     { TrigAddrMask=0xffff, TrigEvent=16, TrigRollover=17 };
enum XilinxTrigCtrl
     { TrigPosMask=0xffff, TrigEnable=16 };
enum XilinxTrigPat
     { TrigPatWord=0xffff, TrigPatShift=16, TrigPort=0 };
enum PortModeFlags
     { BurstMode=0, DataOutMode=1, ClkOutMode=2, ClkSelMode=3,
       ClkModeMask=0x18, IOStartMode=5, HandshakeMode=6,
       IntMode=8, NonDMAMode=9, MemMapped=10,
       DummyIO=12, PortOpen=15, TriggerMode=16,
       PortLargeRam=28, OnePortMode=29, HyperClkMode=30,
       ModeChanged=31
     };
enum PortStatusFlags
     { PortInt=0, PortIntMSB=2, PortIntErr=3,
       PortTrigger=4, PortRollover=5,
       DataAvailable=8, DataIncomplete=9, DataOverrun=10, DataUnderrun=11
     };

typedef struct
{ unsigned short status, block;
  unsigned count;
  long     time_val[2]; /* Equivalent to struct timeval decl. in driver */
  /* double   tv;   32 bits sec, 20 bits usec: exactly rep. in a double */
} BlockType;

enum { KernelBufSize = 1<<17, MaxKernelBufs=32, MaxPortBufs=MaxKernelBufs,
       MaxBlockSize=1<<15, BlockListSize=4096,
       MaxUsrBlocks=BlockListSize/sizeof(BlockType)
     };

typedef struct
{ unsigned flag;
  unsigned short blocks, rblock, wblock, usr_block;
           /* usr_block always refers last rblock in either direction */
  unsigned count, int_overrun, dma_incomplete, overrun, underrun;
  unsigned freq;
  unsigned char *mem_addr, *memp, *mem_end;  /* Xilinx Mem pointers */
  unsigned char *bufp[MaxPortBufs]; /* bufp[i] is a block pointer */
  unsigned mem_size,        /* Xilinx Mem_size */
           buf_ind[MaxPortBufs], /* Default kernel buffer address */
           bufs,            /* Number of kernel bufs for block storage */
           block_size,      /* blocking that is used for a given mode */
                            /* mem_size is always multiple of block_size */
           data_offset;     /* offset within block_size that the user */
                            /* has read or written into */
  unsigned char **blockp;
  BlockType  *block_list;   /* Info about the stored blocks */
} PortInfoType; 

typedef struct { unsigned start, curr, size, flag; } PortAddrType;

typedef struct
{ unsigned short vendor_id, dev_id, command, status;
  unsigned char class[ClassDefs];
  unsigned char cache_line_size, latency_timer, header_type, BIST;
  unsigned base_addr[BaseAddrs];
  unsigned cis_ptr;
  unsigned short sub_vendor_id, sub_sys_id;
  unsigned exp_rom_addr;
  unsigned res[2];
  unsigned char irq_line, irq_pin, min_gnt, max_lat;
} PCIDevConfigType;

int  set_io(unsigned on_off);
int  set_modes(int fd, unsigned enable, unsigned disable);
int  set_port_size(int fd, unsigned size);
int  set_freq(int fd, unsigned clk_mode, unsigned freq);
int  set_trigger(int fd, unsigned pat, unsigned pos);
int  set_one_port(int fd, unsigned on_off);
int  set_extra_bufs(unsigned bufs);
int  set_mmap(int fd, unsigned mode, unsigned size, void **addr);

int  get_pci_conf(int fd, PCIDevConfigType *conf);
int  get_plx_conf(int fd, unsigned *conf, int size);
int  get_xilinx_conf(int fd, unsigned *conf, int size);
int  get_port_info(PortInfoType *portp);
int  get_port_addr(int fd, PortAddrType *addr);
int  get_block_info(BlockType *blockp);

int  read_mem_data(int fd, unsigned words, unsigned offset, unsigned *iobuf);
int  write_mem_data(int fd, unsigned words, unsigned offset, unsigned *iobuf);

void print_status(unsigned short status, FILE *f);
void print_block_info(int fd, FILE *f);
void print_port_info(int fd, FILE *f);
void print_xilinx_regs(int fd, FILE *f);
void print_plx_regs(int fd, FILE *f);
void print_pci_config(int fd, FILE *f);

# endif
