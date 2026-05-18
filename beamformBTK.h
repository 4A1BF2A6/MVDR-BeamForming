#ifndef __BEAMFORMING_BTK_H__
#define __BEAMFORMING_BTK_H__

#define BTK_D (64)

void btk_beamforming_init();
int btk_beamforming_process(short input_channels[][BTK_D],short *out);
void btk_beamforming_set_location(short theta,short phi);
void btk_beamforming_set_location_num(int doa_num);
//void btk_beamforming_set_number(int index);
//int btk_beamforming_process_doa(short input_channels[][64],short *out,int theta_index,int phi_index);

#endif