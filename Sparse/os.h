#ifndef os_h
#define os_h

#if defined _WIN32 || defined __WIN32__ || defined __WINDOWS__ || defined __WIN64__ || defined _WIN64
#define WINDOWS_OS
#elif defined __unix__ || defined __unix
#define UNIXLIKE_OS
#endif

#endif
