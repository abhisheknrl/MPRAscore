/***********************************************************************************
 *
 *   stdout.h -- Standard Output using ITextOutput
 *
 */

#ifndef SYSTEM_STDOUT_H
#define SYSTEM_STDOUT_H

class ITextOutput;
ITextOutput *SetStdout(ITextOutput *pOut);

void system_printf(const char *ach, ...);
void system_putchar(int c);
void system_puts(const char *ach);

// Override by force (may issue warnings with some compilers)
#define printf  system_printf
#define putchar system_putchar
#define puts    system_puts

#endif
