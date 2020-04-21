#!/usr/bin/env python

try:
    from unittest2 import TestCase, main, expectedFailure, skip
except:
    from unittest import TestCase, main, expectedFailure, skip

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG) 
  
from astropy.coordinates import Angle
from datetime import datetime
import trusat.iod as iod

def str_to_bool(s):
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
         raise ValueError # evil ValueError that doesn't tell you what the wrong value was


def test_obs_file_with_answers(obs_type,obs_file,obs_file_answers):
    test_result = True

    with open(obs_file_answers) as file:
        line_number = 0
        answer_key = {}
        last_header = "Beginning of file"
        for line in file:
            line_number += 1
            (answer, test_line) = line.split(',')[0:2]
            if (test_line.startswith("#")):
                last_header = test_line.strip()
            answer_key.update( { line_number : {"answer" : str_to_bool(answer), "last_header" : last_header} } )
            # print(line_number,answer_key[line_number]["answer"],answer_key[line_number]["last_header"])

    with open(obs_file) as file:
        print("Parsing {}".format(obs_file))
        line_number = 0
        for line in file:
            line_number += 1

            if (obs_type is "IOD"):
                result = iod.format_test_iod(line)
            elif (obs_type is "UK"):
                result = iod.format_test_uk(line)
            elif (obs_type is "RDE"):
                result = iod.format_test_rde(line)

            if (result is answer_key[line_number]["answer"]):
                answer = "PASS"
            else:
                answer = "FAIL"
                test_result = False
                print("Line {:>3d} {}: {}\n>>>{}\n".format(
                    line_number, 
                    answer, 
                    answer_key[line_number]["last_header"], 
                    line.rstrip() ))

    return test_result


class Tests(TestCase):
    def test_launch_piece_number_to_letter(self):
        print("launch_piece_number_to_letter(LaunchNum)...")
        A   = iod.launch_piece_number_to_letter('1')
        AA  = iod.launch_piece_number_to_letter('25')
        AAA = iod.launch_piece_number_to_letter('601')
        ZZZ = iod.launch_piece_number_to_letter('14424')

        self.assertEqual(A,'A')
        self.assertEqual(AA,'AA')
        self.assertEqual(AAA,'AAA')
        self.assertEqual(ZZZ,'ZZZ')

    def test_Angle(self):
        print("class Angle...")
        # Format 1: RA/DEC = HHMMSSs+DDMMSS MX   (MX in seconds of arc)
        # 		 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
        # 		 3: RA/DEC = HHMMmmm+DDdddd MX   (MX in degrees of arc)
        # 		 4: AZ/EL  = DDDMMSS+DDMMSS MX   (MX in seconds of arc)
        # 		 5: AZ/EL  = DDDMMmm+DDMMmm MX   (MX in minutes of arc)
        # 		 6: AZ/EL  = DDDdddd+DDdddd MX   (MX in degrees of arc)
        # 		 7: RA/DEC = HHMMSSs+DDdddd MX   (MX in degrees of arc)
        #              IOD 1 0047196+603800
        #              IOD 2 1256405+515478
        #              IOD 3 0859511+119509
        # IOD 4 (20 obs and they look bogus)
        #              IOD 5 1835170-081148 (1855 epoch)
        #              IOD 6 0214030+732900
        #              IOD 7 2047449+293762

        # UK
        #   Code 1: RA/DEC = HHMMSSss+DDMMSSs
        # 	     2: RA/DEC = HHMMmmmm+DDMMmmm
        # 	     3: RA/DEC = HHMMmmmm+DDddddd
        #               UK 1 11203617  +235641 (1950 EPOCH)
        #               UK 2 22202604  +56358
        #               UK 3 13221893  +3777
        # UK 4 (No data)
        # UK 5 (No data)
        # UK 6 (No data)
        # UK 7 (No data)
            # UK 8 (No data)
            # UK 9 (No data) 

        # RDE (Like UK format 1 without subseconds)
        #   Code 1: RA/DEC = HHMMSS+DDMMSS
        #              RDE 1 204836+280552
        #              RDE 1 102009+654911 (1950 EPOCH)


        # Format 1: RA/DEC = HHMMSSs+DDMMSS MX   (MX in seconds of arc)
        #              IOD 1 0047196+603829
        ra_iod1  = Angle('00h47m19.6s')  # (HHMMSSs)
        dec_iod1 = Angle('60d38m29s')    # (DDMMSS)

        angle = iod.Angle(1,5,'0047196+603829',0,"","IOD")
        self.assertAlmostEqual(angle.RA,ra_iod1.deg,7,msg="IOD Format 1 failed in RA")
        self.assertAlmostEqual(angle.DEC,dec_iod1.deg,7, msg="IOD Format 1 failed in DEC")
        self.assertEqual(angle.Epoch,2000)

        # 		 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
        #              IOD 2 1256405+515478
        ra_iod2  = Angle('12h56.405m')  # (HHMMmmm)
        dec_iod2 = Angle('51d54.78m')   # (DDMMmm)

        angle = iod.Angle(2,5,'1256405+515478',0,"","IOD")
        self.assertAlmostEqual(angle.RA,ra_iod2.deg,6, msg="IOD Format 2 failed in RA")
        self.assertAlmostEqual(angle.DEC,dec_iod2.deg,6, msg="IOD Format 2 failed in DEC")
        self.assertEqual(angle.Epoch,2000)

        # 		 3: RA/DEC = HHMMmmm+DDdddd MX   (MX in degrees of arc)
        #              IOD 3 0859511+119509

        ra_iod3  = Angle('08h59.511m')  # (HHMMmmm)
        dec_iod3 = Angle('11.9509d')    # (DDdddd)

        angle = iod.Angle(3,5,'0859511+119509',0,"","IOD")
        self.assertAlmostEqual(angle.RA,ra_iod3.deg,6, msg="IOD Format 3 failed in RA")
        self.assertAlmostEqual(angle.DEC,dec_iod3.deg,6, msg="IOD Format 3 failed in DEC")
        self.assertEqual(angle.Epoch,2000)

        # 		 4: AZ/EL  = DDDMMSS+DDMMSS MX   (MX in seconds of arc)
        #              IOD 4 1835301+603829
        az_iod4  = Angle('183d53m01s')  # (DDDMMSS)
        el_iod4 = Angle('60d38m29s')   # (DDMMSS)

        angle = iod.Angle(4,5,'1835301+603829',0,"","IOD")
        self.assertAlmostEqual(angle.AZ,az_iod4.deg,6, msg="IOD Format 4 failed in AZ")
        self.assertAlmostEqual(angle.EL,el_iod4.deg,6, msg="IOD Format 4 failed in EL")
        self.assertEqual(angle.Epoch,2000)

        # 		 5: AZ/EL  = DDDMMmm+DDMMmm MX   (MX in minutes of arc)
        #              IOD 5 1835170-081148 (1855 epoch)
        az_iod5  = Angle('183d51.70m')  # (DDDMMmm)
        el_iod5 = Angle('-08d11.48m') # (DDMMmm)

        angle = iod.Angle(5,5,'1835170-081148',0,"","IOD")
        self.assertAlmostEqual(angle.AZ,az_iod5.deg,6, msg="IOD Format 5 failed in AZ")
        self.assertAlmostEqual(angle.EL,el_iod5.deg,6, msg="IOD Format 5 failed in EL")
        self.assertEqual(angle.Epoch,2000)

        # 		 6: AZ/EL  = DDDdddd+DDdddd MX   (MX in degrees of arc)
        #              IOD 6 0214030+732900
        az_iod6  = Angle('021.4030d')  # (DDDdddd)
        el_iod6 = Angle('73.2900d')   # (DDdddd)

        angle = iod.Angle(6,5,'0214030+732900',0,"","IOD")
        self.assertAlmostEqual(angle.AZ,az_iod6.deg,6, msg="IOD Format 6 failed in AZ")
        self.assertAlmostEqual(angle.EL,el_iod6.deg,6, msg="IOD Format 6 failed in EL")
        self.assertEqual(angle.Epoch,2000)

        # 		 7: RA/DEC = HHMMSSs+DDdddd MX   (MX in degrees of arc)
        #              IOD 7 2047449+293762
        ra_iod7  = Angle('20h47m44.9s')  # (HHMMSSs)
        dec_iod7 = Angle('29.3762d')     # (DDdddd)

        angle = iod.Angle(7,5,'2047449+293762',0,"","IOD")
        self.assertAlmostEqual(angle.RA,ra_iod7.deg,6, msg="IOD Format 7 failed in RA")
        self.assertAlmostEqual(angle.DEC,dec_iod7.deg,6, msg="IOD Format 7 failed in DEC")
        self.assertEqual(angle.Epoch,2000)

        # RDE
        #   Code 1: RA/DEC = HHMMSSss+DDMMSSs
        #              RDE 1 202140+005931

        ra_rde1  = Angle('20h21m40s')  # (HHMMSSss)
        dec_rde1 = Angle('00d59m31s')  # (DDMMSSs)

        angle = iod.Angle(1,5,'202140+005931',0,"","RDE")
        self.assertAlmostEqual(angle.RA,ra_rde1.deg,6, msg="RDE Format 1 failed in RA")
        self.assertAlmostEqual(angle.DEC,dec_rde1.deg,6, msg="RDE Format 1 failed in DEC")
        self.assertEqual(angle.Epoch,2000)

        #UK
        #               UK 1 11203617  +235641 (1950 EPOCH)

        ra_uk1  = Angle('11h20m36.17s')  # (HHMMSSss)
        dec_uk1 = Angle('23d56m41s')     # (DDMMSSs)

        angle = iod.Angle(1,5,'11203617+235641',0,"","UK")
        self.assertAlmostEqual(angle.RA,ra_uk1.deg,6, msg="UK Format 1 failed in RA")
        self.assertAlmostEqual(angle.DEC,dec_uk1.deg,6, msg="UK Format 1 failed in DEC")
        self.assertEqual(angle.Epoch,2000)

        # 	     2: RA/DEC = HHMMmmmm+DDMMmmm
        #               UK 2 22202604  +56358
        ra_uk2  = Angle('22h20.2604m')  # (HHMMmmmm)
        dec_uk2 = Angle('56d35.8m')     # (DDMMmmm)

        angle = iod.Angle(2,5,'22202604+56358',0,"","UK")
        self.assertAlmostEqual(angle.RA,ra_uk2.deg,6, msg="UK Format 2 failed in RA")
        self.assertAlmostEqual(angle.DEC,dec_uk2.deg,6, msg="UK Format 2 failed in DEC")
        self.assertEqual(angle.Epoch,2000)

        # 	     3: RA/DEC = HHMMmmmm+DDddddd
        #               UK 3 13221893  +3777
        ra_uk3  = Angle('13h22.1893m')  # (HHMMmmmm)
        dec_uk3 = Angle('37.77d')       # (DDddddd)

        angle = iod.Angle(3,5,'13221893+3777',0,"","UK")
        self.assertAlmostEqual(angle.RA,ra_uk3.deg,6, msg="UK Format 3 failed in RA")
        self.assertAlmostEqual(angle.DEC,dec_uk3.deg,6, msg="UK Format 3 failed in DEC")
        self.assertEqual(angle.Epoch,2000)


        ###
        ### Tests that drop significant figures from the right
        ###
        # HHMMSSs
        HHMMSS_ts = iod.angle_from_HHMMSSss('0047191')
        HHMMSS_ap = Angle('00h47m19.1s')   # (HHMMSSs)
        self.assertAlmostEqual(HHMMSS_ts,HHMMSS_ap.deg,7,msg="7 SIGFIG HHMMSSs ")

        HHMMSS_ts = iod.angle_from_HHMMSSss('004719')
        HHMMSS_ap = Angle('00h47m19s')   # (HHMMSSs)
        self.assertAlmostEqual(HHMMSS_ts,HHMMSS_ap.deg,7,msg="6 SIGFIG HHMMSSs ")

        HHMMSS_ts = iod.angle_from_HHMMSSss('00471')
        HHMMSS_ap = Angle('00h47m10s')   # (HHMMSSs)
        self.assertAlmostEqual(HHMMSS_ts,HHMMSS_ap.deg,7,msg="5 SIGFIG HHMMSSs ")

        HHMMSS_ts = iod.angle_from_HHMMSSss('0047')
        HHMMSS_ap = Angle('00h47m')   # (HHMMSSs)
        self.assertAlmostEqual(HHMMSS_ts,HHMMSS_ap.deg,7,msg="4 SIGFIG HHMMSSs ")

        HHMMSS_ts = iod.angle_from_HHMMSSss('004')
        HHMMSS_ap = Angle('00h40m')   # (HHMMSSs)
        self.assertAlmostEqual(HHMMSS_ts,HHMMSS_ap.deg,7,msg="3 SIGFIG HHMMSSs ")

        HHMMSS_ts = iod.angle_from_HHMMSSss('09')
        HHMMSS_ap = Angle('09h')   # (HHMMSSs)
        self.assertAlmostEqual(HHMMSS_ts,HHMMSS_ap.deg,7,msg="2 SIGFIG HHMMSSs ")

        HHMMSS_ts = iod.angle_from_HHMMSSss('1')
        HHMMSS_ap = Angle('10h')   # (HHMMSSs)
        self.assertAlmostEqual(HHMMSS_ts,HHMMSS_ap.deg,7,msg="1 SIGFIG HHMMSSs ")

        # 13221893 HHMMmmmm

        HHMMmmmm_ts = iod.angle_from_HHMMmmmm('13221893')
        HHMMmmmm_ap = Angle('13h22.1893m')  # (HHMMmmmm)
        self.assertAlmostEqual(HHMMmmmm_ts,HHMMmmmm_ap.deg,7,msg="8 SIGFIG HHMMmmmm")

        HHMMmmmm_ts = iod.angle_from_HHMMmmmm('1322189')
        HHMMmmmm_ap = Angle('13h22.189m')  # (HHMMmmmm)
        self.assertAlmostEqual(HHMMmmmm_ts,HHMMmmmm_ap.deg,7,msg="7 SIGFIG HHMMmmmm")

        HHMMmmmm_ts = iod.angle_from_HHMMmmmm('132218')
        HHMMmmmm_ap = Angle('13h22.18m')  # (HHMMmmmm)
        self.assertAlmostEqual(HHMMmmmm_ts,HHMMmmmm_ap.deg,7,msg="6 SIGFIG HHMMmmmm")

        HHMMmmmm_ts = iod.angle_from_HHMMmmmm('13221')
        HHMMmmmm_ap = Angle('13h22.1m')  # (HHMMmmmm)
        self.assertAlmostEqual(HHMMmmmm_ts,HHMMmmmm_ap.deg,7,msg="5 SIGFIG HHMMmmmm")

        HHMMmmmm_ts = iod.angle_from_HHMMmmmm('1322')
        HHMMmmmm_ap = Angle('13h22m')  # (HHMMmmmm)
        self.assertAlmostEqual(HHMMmmmm_ts,HHMMmmmm_ap.deg,7,msg="4 SIGFIG HHMMmmmm")

        HHMMmmmm_ts = iod.angle_from_HHMMmmmm('132')
        HHMMmmmm_ap = Angle('13h20m')  # (HHMMmmmm)
        self.assertAlmostEqual(HHMMmmmm_ts,HHMMmmmm_ap.deg,7,msg="3 SIGFIG HHMMmmmm")

        HHMMmmmm_ts = iod.angle_from_HHMMmmmm('13')
        HHMMmmmm_ap = Angle('13h')  # (HHMMmmmm)
        self.assertAlmostEqual(HHMMmmmm_ts,HHMMmmmm_ap.deg,7,msg="2 SIGFIG HHMMmmmm")

        HHMMmmmm_ts = iod.angle_from_HHMMmmmm('1')
        HHMMmmmm_ap = Angle('10h')  # (HHMMmmmm)
        self.assertAlmostEqual(HHMMmmmm_ts,HHMMmmmm_ap.deg,7,msg="1 SIGFIG HHMMmmmm")

        # 18353513 # (DDDMMSSs)

        DDDMMSSs_ts = iod.angle_from_DDDMMSSs('18353513')
        DDDMMSSs_ap = Angle('183d53m51.3s')  # (DDDMMSSs)
        self.assertAlmostEqual(DDDMMSSs_ts,DDDMMSSs_ap.deg,7,msg="8 SIGFIG DDDMMSSs")

        DDDMMSSs_ts = iod.angle_from_DDDMMSSs('1835351')
        DDDMMSSs_ap = Angle('183d53m51s')  # (DDDMMSSs)
        self.assertAlmostEqual(DDDMMSSs_ts,DDDMMSSs_ap.deg,7,msg="7 SIGFIG DDDMMSSs")

        DDDMMSSs_ts = iod.angle_from_DDDMMSSs('183535')
        DDDMMSSs_ap = Angle('183d53m50s')  # (DDDMMSSs)
        self.assertAlmostEqual(DDDMMSSs_ts,DDDMMSSs_ap.deg,7,msg="6 SIGFIG DDDMMSSs")

        DDDMMSSs_ts = iod.angle_from_DDDMMSSs('18353')
        DDDMMSSs_ap = Angle('183d53m')  # (DDDMMSSs)
        self.assertAlmostEqual(DDDMMSSs_ts,DDDMMSSs_ap.deg,7,msg="5 SIGFIG DDDMMSSs")

        DDDMMSSs_ts = iod.angle_from_DDDMMSSs('1835')
        DDDMMSSs_ap = Angle('183d50m')  # (DDDMMSSs)
        self.assertAlmostEqual(DDDMMSSs_ts,DDDMMSSs_ap.deg,7,msg="4 SIGFIG DDDMMSSs")

        DDDMMSSs_ts = iod.angle_from_DDDMMSSs('183')
        DDDMMSSs_ap = Angle('183d')  # (DDDMMSSs)
        self.assertAlmostEqual(DDDMMSSs_ts,DDDMMSSs_ap.deg,7,msg="3 SIGFIG DDDMMSSs")

        DDDMMSSs_ts = iod.angle_from_DDDMMSSs('18')
        DDDMMSSs_ap = Angle('180d')  # (DDDMMSSs)
        self.assertAlmostEqual(DDDMMSSs_ts,DDDMMSSs_ap.deg,7,msg="2 SIGFIG DDDMMSSs")

        DDDMMSSs_ts = iod.angle_from_DDDMMSSs('1')
        DDDMMSSs_ap = Angle('100d')  # (DDDMMSSs)
        self.assertAlmostEqual(DDDMMSSs_ts,DDDMMSSs_ap.deg,7,msg="1 SIGFIG DDDMMSSs")

        # +156358654 # (DDDMMmmm)

        DDDMMmmm_ts = iod.angle_from_DDDMMmmm('15635865')
        DDDMMmmm_ap = Angle('156d35.865m')  # (DDDMMmmm)
        self.assertAlmostEqual(DDDMMmmm_ts,DDDMMmmm_ap.deg,7,msg="8 SIGFIG DDDMMmmm")

        DDDMMmmm_ts = iod.angle_from_DDDMMmmm('1563586')
        DDDMMmmm_ap = Angle('156d35.86m')  # (DDDMMmmm)
        self.assertAlmostEqual(DDDMMmmm_ts,DDDMMmmm_ap.deg,7,msg="7 SIGFIG DDDMMmmm")

        DDDMMmmm_ts = iod.angle_from_DDDMMmmm('156358')
        DDDMMmmm_ap = Angle('156d35.8m')  # (DDDMMmmm)
        self.assertAlmostEqual(DDDMMmmm_ts,DDDMMmmm_ap.deg,7,msg="6 SIGFIG DDDMMmmm")

        DDDMMmmm_ts = iod.angle_from_DDDMMmmm('15635')
        DDDMMmmm_ap = Angle('156d35m')  # (DDDMMmmm)
        self.assertAlmostEqual(DDDMMmmm_ts,DDDMMmmm_ap.deg,7,msg="5 SIGFIG DDDMMmmm")

        DDDMMmmm_ts = iod.angle_from_DDDMMmmm('1563')
        DDDMMmmm_ap = Angle('156d30m')  # (DDDMMmmm)
        self.assertAlmostEqual(DDDMMmmm_ts,DDDMMmmm_ap.deg,7,msg="4 SIGFIG DDDMMmmm")

        DDDMMmmm_ts = iod.angle_from_DDDMMmmm('156')
        DDDMMmmm_ap = Angle('156d')  # (DDDMMmmm)
        self.assertAlmostEqual(DDDMMmmm_ts,DDDMMmmm_ap.deg,7,msg="3 SIGFIG DDDMMmmm")

        DDDMMmmm_ts = iod.angle_from_DDDMMmmm('15')
        DDDMMmmm_ap = Angle('150d')  # (DDDMMmmm)
        self.assertAlmostEqual(DDDMMmmm_ts,DDDMMmmm_ap.deg,7,msg="2 SIGFIG DDDMMmmm")

        DDDMMmmm_ts = iod.angle_from_DDDMMmmm('1')
        DDDMMmmm_ap = Angle('100d')  # (DDDMMmmm)
        self.assertAlmostEqual(DDDMMmmm_ts,DDDMMmmm_ap.deg,7,msg="1 SIGFIG DDDMMmmm")

        # 17854321 # (DDDddddd)
        DDDddddd_ts = iod.angle_from_DDDddddd('17854321')
        DDDddddd_ap = Angle('178.54321d')  # ('178.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="8 SIGFIG DDDddddd")

        DDDddddd_ts = iod.angle_from_DDDddddd('1785432')
        DDDddddd_ap = Angle('178.5432d')  # ('178.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="7 SIGFIG DDDddddd")

        DDDddddd_ts = iod.angle_from_DDDddddd('178543')
        DDDddddd_ap = Angle('178.543d')  # ('178.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="6 SIGFIG DDDddddd")

        DDDddddd_ts = iod.angle_from_DDDddddd('178543')
        DDDddddd_ap = Angle('178.543d')  # ('178.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="6 SIGFIG DDDddddd")

        DDDddddd_ts = iod.angle_from_DDDddddd('17854')
        DDDddddd_ap = Angle('178.54d')  # ('178.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="5 SIGFIG DDDddddd")

        DDDddddd_ts = iod.angle_from_DDDddddd('1785')
        DDDddddd_ap = Angle('178.5d')  # ('178.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="4 SIGFIG DDDddddd")

        DDDddddd_ts = iod.angle_from_DDDddddd('178')
        DDDddddd_ap = Angle('178d')  # ('178.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="3 SIGFIG DDDddddd")

        DDDddddd_ts = iod.angle_from_DDDddddd('17')
        DDDddddd_ap = Angle('170d')  # ('178.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="2 SIGFIG DDDddddd")

        DDDddddd_ts = iod.angle_from_DDDddddd('1')
        DDDddddd_ap = Angle('100d')  # ('178.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="1 SIGFIG DDDddddd")

        # Negative angles

        DDDMMSSs_ts = iod.angle_from_DDDMMSSs('-835351')
        DDDMMSSs_ap = Angle('-83d53m51s')  # (DDDMMSSs)
        self.assertAlmostEqual(DDDMMSSs_ts,DDDMMSSs_ap.deg,7,msg="7 SIGFIG -DDMMSSs")

        DDDMMmmm_ts = iod.angle_from_DDDMMmmm('-563586')
        DDDMMmmm_ap = Angle('-56d35.86m')  # (-DDMMmmm)
        self.assertAlmostEqual(DDDMMmmm_ts,DDDMMmmm_ap.deg,7,msg="7 SIGFIG -DDMMmmm")

        DDDddddd_ts = iod.angle_from_DDDddddd('-785432')
        DDDddddd_ap = Angle('-78.5432d')  # ('-78.54321')
        self.assertAlmostEqual(DDDddddd_ts,DDDddddd_ap.deg,7,msg="7 SIGFIG -DDddddd")

    def test_get_angle_content_IOD(self):
        print("get_angle_content_IOD(Anglestring)...")
        angle_cases = [
            (0,  '0047196+603829',  True),  # Standard format
            (1,  '0859511+1199',    True),  # Fewer sig figs in 1st component
            (2,  '18353  -603829',  True),  # Trailing space in 1st component
            (3,  '0047196+6038  ',  True),  # Trailing space in 2nd component
            (4,  '183    -081   ',  True),  # Trailing space in both components
            (5,    '12564+515478',  False), # Fewer sig figs in 1st component
            (6,  ' 047449+293762',  False), # Leading space in 1st component
            (7,  'A047196+603829',  False), # Character in 1st component
            (8,  '1256405 515478',  False), # Missing sign
            (9,  '1047196+ 63829',  False), # Leading space in 2nd component
            (10, '2047449+293.76',  False), # Decimal point in 2nd component
            (11, '00471962+603829', False), # Too many digits in 1st component
            (12, '0047196+6038291', False), # Too many digits in 2nd component
            (13, '0047 96+603829',  False), # Space in middle of 1st component
            (14, '0047196+6038 9',  False)  # Space in middle of 2nd component
        ]

        for (case_no, AngleString, answer) in angle_cases:
            (Angle1,Angle2) = iod.get_angle_content_IOD(AngleString)
            if (Angle1 and Angle1):
                result = True
            else:
                result = False
            self.assertEqual(result,answer,msg="IOD Case {}  AngleString='{}'".format(case_no,AngleString))

    def test_get_angle_content_UK(self):
        print("get_angle_content_UK(Anglestring)...")
        angle_cases = [
            (0,  '00471961+6038291',  True),  # Standard format
            (1,  '08595111+1199',     True),  # Fewer sig figs in 1st component
            (2,  '183532  -6038291',  True),  # Trailing space in 1st component
            (3,  '00471961+6038   ',  True),  # Trailing space in 2nd component
            (4,  '183     -081    ',  True),  # Trailing space in both components
            (5,    '125641+5154781',  False), # Fewer sig figs in 1st component
            (6,  ' 0474491+2937621',  False), # Leading space in 1st component
            (7,  'A0471961+6038291',  False), # Character in 1st component
            (8,  '12564051 5154781',  False), # Missing sign
            (9,  '10471961+ 638291',  False), # Leading space in 2nd component
            (10, '20474491+293.761',  False), # Decimal point in 2nd component
            (11, '004719621+6038291', False), # Too many digits in 1st component
            (12, '00471961+60382911', False), # Too many digits in 2nd component
            (13, '0047 961+6038291',  False), # Space in middle of 1st component
            (14, '00471961+6038 91',  False)  # Space in middle of 2nd component
        ]

        for (case_no, AngleString, answer) in angle_cases:
            (Angle1,Angle2) = iod.get_angle_content_UK(AngleString)
            if (Angle1 and Angle1):
                result = True
            else:
                result = False
            self.assertEqual(result,answer,msg="UK Case {}  AngleString='{}'".format(case_no,AngleString))

    def test_get_angle_content_RDE(self):
        print("get_angle_content_RDE(Anglestring)...")

        angle_cases = [
            (0,  '004719+603829',  True),  # Standard format
            (1,  '085951+1199',    True),  # Fewer sig figs  in 1st component
            (2,  '18353 -603829',  True),  # Trailing space  in 1st component
            (3,  '004719+6038  ',  True),  # Trailing space  in 2nd component
            (4,  '183   -081   ',  True),  # Trailing space  in both components
            (5,    '1256+515478',  False), # Fewer sig figs  in 1st component
            (6,  ' 04744+293762',  False), # Leading space   in 1st component
            (7,  'A04719+603829',  False), # Character       in 1st component
            (8,  '125640 515478',  False), # Missing sign
            (9,  '104719+ 63829',  False), # Leading space   in 2nd component
            (10, '204744+293.76',  False), # Decimal point   in 2nd component
            (11, '0047191+603829', False), # Too many digits in 1st component
            (12, '004719+6038291', False), # Too many digits in 2nd component
            (13, '0047 9+603829',  False), # Space in middle of 1st component
            (14, '004719+6038 9',  False)  # Space in middle of 2nd component
        ]

        for (case_no, AngleString, answer) in angle_cases:
            (Angle1,Angle2) = iod.get_angle_content_RDE(AngleString)
            if (Angle1 and Angle1):
                result = True
            else:
                result = False
            self.assertEqual(result,answer,msg="RDE Case {}  AngleString='{}'".format(case_no,AngleString))

    def test_format_test_iod_compliant(self):
        print("format_test_iod(block)...")

        obs_file = "data/tests/IOD_data_compliant.txt" 
        obs_file_answers = "data/tests/IOD_data_compliant_answers.txt"
        obs_type = "IOD"
        self.assertEqual(test_obs_file_with_answers(obs_type,obs_file,obs_file_answers),True,msg="IOD compliant test failed.")


    @expectedFailure
    def test_format_test_iod_noncompliant(self):
        print("format_test_iod(non compliant block)...")

        obs_file = "data/tests/IOD_data_noncompliant.txt" 
        obs_file_answers = "data/tests/IOD_data_noncompliant_answers.txt"
        obs_type = "IOD"
        self.assertEqual(test_obs_file_with_answers(obs_type,obs_file,obs_file_answers),True,msg="IOD non compliant test failed.")


    def test_parse_iod_lines(self):
        print("test_parse_iod_lines(block)...")
        test_file = "data/tests/IOD_data_compliant.txt" 
        IODs = iod.get_iod_records_from_file(test_file)

        self.assertEqual(len(IODs),46,msg="Failed to read all (46) compliant IODs")

        test_iod_string = "23794 96 010A   2701 G 20040506012614270 17 25 1100114-184298 38 I+020 10  10000 # Comments after the line here"
        IOD = iod.parse_iod_lines(test_iod_string)[0]
        self.assertEqual(IOD.ObjectNumber, 23794 ,msg="Failed ObjectNumber {}".format(IOD.ObjectNumber))
        self.assertEqual(IOD.LaunchYear, 1996 ,msg="Failed LaunchYear {}".format(IOD.LaunchYear))
        self.assertEqual(IOD.InternationalDesignation, "1996-010A  " ,msg="Failed InternationalDesignation {}".format(IOD.InternationalDesignation))
        self.assertEqual(IOD.Station, 2701 ,msg="Failed Station {}".format(IOD.Station))
        self.assertEqual(IOD.StationStatusCode, "G" ,msg="Failed StationStatusCode {}".format(IOD.StationStatusCode))
        self.assertEqual(IOD.StationStatus, None ,msg="Failed StationStatus {}".format(IOD.StationStatus))
        self.assertEqual(IOD.DateTimeString, "20040506012614270" ,msg="Failed DateTimeString {}".format(IOD.DateTimeString))
        self.assertEqual(IOD.DateTime, datetime(2004,5,6,1,26,14,270000) ,msg="Failed DateTime {}".format(IOD.DateTime))
        self.assertEqual(IOD.TimeUncertainty,0.1 ,msg="Failed TimeUncertainty {}".format(IOD.TimeUncertainty))
        self.assertEqual(IOD.TimeStandardCode,None ,msg="Failed TimeStandardCode {}".format(IOD.TimeStandardCode))
        self.assertEqual(IOD.AngleFormatCode,2 ,msg="Failed AngleFormatCode {}".format(IOD.AngleFormatCode))
        self.assertEqual(IOD.EpochCode,5 ,msg="Failed EpochCode {}".format(IOD.EpochCode))
        self.assertEqual(IOD.Epoch,2000 ,msg="Failed Epoch {}".format(IOD.Epoch))
        self.assertEqual(IOD.RA,165.02849999999998,msg="Failed RA {}".format(IOD.RA))
        self.assertEqual(IOD.DEC,-18.71633333333333 ,msg="Failed DEC {}".format(IOD.DEC))
        self.assertEqual(IOD.AZ,0 ,msg="Failed AZ {}".format(IOD.AZ))
        self.assertEqual(IOD.EL,0 ,msg="Failed EL {}".format(IOD.EL))
        self.assertEqual(IOD.ValidPosition,1 ,msg="Failed ValidPosition {}".format(IOD.ValidPosition))
        self.assertEqual(IOD.PositionUncertainty,0.05 ,msg="Failed PositionUncertainty {}".format(IOD.PositionUncertainty))
        self.assertEqual(IOD.OpticalCode,"I",msg="Failed OpticalCode {}".format(IOD.OpticalCode))
        self.assertEqual(IOD.VisualMagnitude,2.0 ,msg="Failed VisualMagnitude {}".format(IOD.VisualMagnitude))
        self.assertEqual(IOD.MagnitudeUncertainty,1.0 ,msg="Failed MagnitudeUncertainty {}".format(IOD.MagnitudeUncertainty))
        self.assertEqual(IOD.FlashPeriod,10.0 ,msg="Failed FlashPeriod {}".format(IOD.FlashPeriod))
        self.assertEqual(IOD.VisualMagnitude_high,None ,msg="Failed VisualMagnitude_high {}".format(IOD.VisualMagnitude_high))
        self.assertEqual(IOD.VisualMagnitude_low,None ,msg="Failed VisualMagnitude_low {}".format(IOD.VisualMagnitude_low))
        self.assertEqual(IOD.Remarks," # Comments after the line here" ,msg="Failed Remarks {}".format(IOD.Remarks))
        self.assertEqual(IOD.message_id,None ,msg="Failed message_id {}".format(IOD.message_id))
        self.assertEqual(IOD.IODType,"IOD" ,msg="Failed IODType {}".format(IOD.IODType))
        self.assertEqual(IOD.iod_string,test_iod_string ,msg="Failed iod_string {}".format(IOD.iod_string))


    # Our regex isn't perfect, some badly formatted records get through the first screen
    # @expectedFailure
    def test_parse_iod_lines_noncompliant(self):
        print("test_parse_iod_lines(block)...")

        test_file = "data/tests/IOD_data_noncompliant.txt" 
        IODs = iod.get_iod_records_from_file(test_file)
        print("Parsed {} lines from {}".format(len(IODs),test_file))


    def test_format_test_uk_compliant(self):
        print("format_test_uk(block)...")

        obs_file = "data/tests/UK_data_compliant.txt" 
        obs_file_answers = "data/tests/UK_data_compliant_answers.txt"
        obs_type = "UK"
        self.assertEqual(test_obs_file_with_answers(obs_type,obs_file,obs_file_answers),True,msg="UK compliant test failed.")

    # Our regex isn't perfect, some badly formatted records get through the first screen
    @expectedFailure
    def test_format_test_uk_noncompliant(self):
        print("format_test_uk(non compliant block)...")

        obs_file = "data/tests/UK_data_noncompliant.txt" 
        obs_file_answers = "data/tests/UK_data_noncompliant_answers.txt"
        obs_type = "UK"
        self.assertEqual(test_obs_file_with_answers(obs_type,obs_file,obs_file_answers),True,msg="UK non compliant test failed.")


    def test_parse_uk_lines(self):
        print("test_parse_uk_lines(block)...")

        obs_file = "data/tests/UK_data_compliant.txt" 
        IODs = iod.get_uk_records_from_file(obs_file)
        self.assertEqual(len(IODs),49,msg="Failed to read all (49) compliant IODs")

        test_uk_string = "9701201201803101520195542  01   12172038  +15585   1  5             +6 +8   190R# Comments after the line here"
        # test_uk_string = "9701201201803101520195542  01   12172038  +15585   1  5             +6 +8   190R# Comments after the line here"
        #                   00000000011111111112222222222333333333344444444445555555555666666666677777777778
        #                   12345678901234567890123456789012345678901234567890123456789012345678901234567890

        IOD = iod.parse_uk_lines(test_uk_string)[0]
        
        self.assertEqual(IOD.ObjectNumber, 0 ,msg="Failed ObjectNumber {}".format(IOD.ObjectNumber))
        self.assertEqual(IOD.LaunchYear, 1997 ,msg="Failed LaunchYear {}".format(IOD.LaunchYear))
        self.assertEqual(IOD.InternationalDesignation, "1997-012A" ,msg="Failed InternationalDesignation {}".format(IOD.InternationalDesignation))
        self.assertEqual(IOD.Station, 2018 ,msg="Failed Station {}".format(IOD.Station))
        self.assertEqual(IOD.StationStatusCode, None ,msg="Failed StationStatusCode {}".format(IOD.StationStatusCode))
        self.assertEqual(IOD.StationStatus, None ,msg="Failed StationStatus {}".format(IOD.StationStatus))
        self.assertEqual(IOD.DateTimeString, "03101520195542  " ,msg="Failed DateTimeString {}".format(IOD.DateTimeString))
        self.assertEqual(IOD.DateTime, datetime(2003, 10, 15, 20, 19, 55, 420000) ,msg="Failed DateTime {}".format(IOD.DateTime))
        self.assertEqual(IOD.TimeUncertainty,0.1 ,msg="Failed TimeUncertainty {}".format(IOD.TimeUncertainty))
        self.assertEqual(IOD.TimeStandardCode,1 ,msg="Failed TimeStandardCode {}".format(IOD.TimeStandardCode))
        self.assertEqual(IOD.AngleFormatCode,2 ,msg="Failed AngleFormatCode {}".format(IOD.AngleFormatCode))
        self.assertEqual(IOD.EpochCode,5 ,msg="Failed EpochCode {}".format(IOD.EpochCode))
        self.assertEqual(IOD.Epoch,2000 ,msg="Failed Epoch {}".format(IOD.Epoch))
        self.assertEqual(IOD.RA,260.09499999999997,msg="Failed RA {}".format(IOD.RA))
        self.assertEqual(IOD.DEC,15.975 ,msg="Failed DEC {}".format(IOD.DEC))
        self.assertEqual(IOD.AZ,0 ,msg="Failed AZ {}".format(IOD.AZ))
        self.assertEqual(IOD.EL,0 ,msg="Failed EL {}".format(IOD.EL))
        self.assertEqual(IOD.ValidPosition,1 ,msg="Failed ValidPosition {}".format(IOD.ValidPosition))
        self.assertEqual(IOD.PositionUncertainty,0.016666666666666666 ,msg="Failed PositionUncertainty {}".format(IOD.PositionUncertainty))
        self.assertEqual(IOD.OpticalCode,"R",msg="Failed OpticalCode {}".format(IOD.OpticalCode))
        self.assertEqual(IOD.VisualMagnitude,None ,msg="Failed VisualMagnitude {}".format(IOD.VisualMagnitude))
        self.assertEqual(IOD.MagnitudeUncertainty,None ,msg="Failed MagnitudeUncertainty {}".format(IOD.MagnitudeUncertainty))
        self.assertEqual(IOD.FlashPeriod,1.9 ,msg="Failed FlashPeriod {}".format(IOD.FlashPeriod))
        self.assertEqual(IOD.VisualMagnitude_high,6 ,msg="Failed VisualMagnitude_high {}".format(IOD.VisualMagnitude_high))
        self.assertEqual(IOD.VisualMagnitude_low,8.0 ,msg="Failed VisualMagnitude_low {}".format(IOD.VisualMagnitude_low))
        self.assertEqual(IOD.Remarks,"# Comments after the line here" ,msg="Failed Remarks {}".format(IOD.Remarks))
        self.assertEqual(IOD.message_id,None ,msg="Failed message_id {}".format(IOD.message_id))
        self.assertEqual(IOD.IODType,"UK" ,msg="Failed IODType {}".format(IOD.IODType))
        self.assertEqual(IOD.iod_string,test_uk_string ,msg="Failed iod_string {}".format(IOD.iod_string))


    def test_get_rde_records_from_file(self):
        print("get_rde_records_from_file(file)...")

        obs_file = "data/tests/RDE_data_compliant.txt" 
        IODs = iod.get_rde_records_from_file(obs_file)
        self.assertEqual(len(IODs),14,msg="RDE compliant test failed.")

    @skip("Requires completion of write_uk_line()")
    def test_parse_rde_block(self):

        test_rde_string="2420 1909 0.211 1204\n28\n0502402 184758.93 193923+452358 3.1 3.1 0 S\n999"

        IOD = iod.parse_rde_record(test_rde_string)
        
        # TODO - Format the RDE "line" as a UK-formatted string.  Work ongoing in iod.py
        iod.write_uk_line(IOD)

        self.assertEqual(IOD.ObjectNumber, 0 ,msg="Failed ObjectNumber {}".format(IOD.ObjectNumber))
        self.assertEqual(IOD.LaunchYear, 2005 ,msg="Failed LaunchYear {}".format(IOD.LaunchYear))
        self.assertEqual(IOD.InternationalDesignation, "2005-024B" ,msg="Failed InternationalDesignation {}".format(IOD.InternationalDesignation))
        self.assertEqual(IOD.Station, 2420 ,msg="Failed Station {}".format(IOD.Station))
        self.assertEqual(IOD.StationStatusCode, None ,msg="Failed StationStatusCode {}".format(IOD.StationStatusCode))
        self.assertEqual(IOD.StationStatus, None ,msg="Failed StationStatus {}".format(IOD.StationStatus))
        self.assertEqual(IOD.DateTimeString, "190928184758.93" ,msg="Failed DateTimeString {}".format(IOD.DateTimeString))
        self.assertEqual(IOD.DateTime, datetime(2019, 9, 28, 18, 47, 58, 930000) ,msg="Failed DateTime {}".format(IOD.DateTime))
        self.assertEqual(IOD.TimeUncertainty,0.2 ,msg="Failed TimeUncertainty {}".format(IOD.TimeUncertainty))
        self.assertEqual(IOD.TimeStandardCode,1 ,msg="Failed TimeStandardCode {}".format(IOD.TimeStandardCode))
        self.assertEqual(IOD.AngleFormatCode,1 ,msg="Failed AngleFormatCode {}".format(IOD.AngleFormatCode))
        self.assertEqual(IOD.EpochCode,4 ,msg="Failed EpochCode {}".format(IOD.EpochCode))
        self.assertEqual(IOD.Epoch,1950 ,msg="Failed Epoch {}".format(IOD.Epoch))
        self.assertEqual(IOD.RA,294.8458333333333,msg="Failed RA {}".format(IOD.RA))
        self.assertEqual(IOD.DEC,45.39944444444444 ,msg="Failed DEC {}".format(IOD.DEC))
        self.assertEqual(IOD.AZ,0 ,msg="Failed AZ {}".format(IOD.AZ))
        self.assertEqual(IOD.EL,0 ,msg="Failed EL {}".format(IOD.EL))
        self.assertEqual(IOD.ValidPosition,1 ,msg="Failed ValidPosition {}".format(IOD.ValidPosition))
        self.assertEqual(IOD.PositionUncertainty,0.03333333333333333 ,msg="Failed PositionUncertainty {}".format(IOD.PositionUncertainty))
        self.assertEqual(IOD.OpticalCode,"S",msg="Failed OpticalCode {}".format(IOD.OpticalCode))
        self.assertEqual(IOD.VisualMagnitude,None ,msg="Failed VisualMagnitude {}".format(IOD.VisualMagnitude))
        self.assertEqual(IOD.MagnitudeUncertainty,None ,msg="Failed MagnitudeUncertainty {}".format(IOD.MagnitudeUncertainty))
        self.assertEqual(IOD.FlashPeriod,0 ,msg="Failed FlashPeriod {}".format(IOD.FlashPeriod))
        self.assertEqual(IOD.VisualMagnitude_high,3.1 ,msg="Failed VisualMagnitude_high {}".format(IOD.VisualMagnitude_high))
        self.assertEqual(IOD.VisualMagnitude_low,3.1 ,msg="Failed VisualMagnitude_low {}".format(IOD.VisualMagnitude_low))
        self.assertEqual(IOD.Remarks,"" ,msg="Failed Remarks {}".format(IOD.Remarks))
        self.assertEqual(IOD.message_id,None ,msg="Failed message_id {}".format(IOD.message_id))
        self.assertEqual(IOD.IODType,"RDE" ,msg="Failed IODType {}".format(IOD.IODType))
        # Not currently formatted correctly
        # self.assertEqual(IOD.iod_string,test_rde_string ,msg="Failed iod_string {}".format(IOD.iod_string))        


    def test_DateTime_frompacked(self):
        print("DateTime_frompacked(DateTimeString,format_type)...")
        # All should test True
        time_cases = [
            (0,   '9909101234567891', 'UK', datetime(1999,9,10,12,34,56,789100) ),
            (1,   '990910123456789 ', 'UK', datetime(1999,9,10,12,34,56,789000) ),
            (2,   '99091012345678  ', 'UK', datetime(1999,9,10,12,34,56,780000) ),
            (3,   '9909101234567   ', 'UK', datetime(1999,9,10,12,34,56,700000) ),
            (4,   '990910123456    ', 'UK', datetime(1999,9,10,12,34,56,000000) ),
            (5,   '9909101234      ', 'UK', datetime(1999,9,10,12,34,00,000000) ),
            (6,   '99091012        ', 'UK', datetime(1999,9,10,12,00,00,000000) ),
            (7,   '990910          ', 'UK', datetime(1999,9,10,00,00,00,000000) ),
            (8,   '990910',           'UK', datetime(1999,9,10,00,00,00,000000) ),
            (9, '19990910123456789', 'IOD', datetime(1999,9,10,12,34,56,789000) ),
            (10,'1999091012345678 ', 'IOD', datetime(1999,9,10,12,34,56,780000) ),
            (11,'199909101234567  ', 'IOD', datetime(1999,9,10,12,34,56,700000) ),
            (12,'19990910123456   ', 'IOD', datetime(1999,9,10,12,34,56,000000) ),
            (13,'199909101234     ', 'IOD', datetime(1999,9,10,12,34,00,000000) ),
            (14,'1999091012       ', 'IOD', datetime(1999,9,10,12,00,00,000000) ),
            (15,'19990910         ', 'IOD', datetime(1999,9,10,00,00,00,000000) ),
            (16,'19990910',          'IOD', datetime(1999,9,10,00,00,00,000000) ),
            (17,  '990910123456.78', 'RDE', datetime(1999,9,10,12,34,56,780000) ),
            (18,  '990910123456.7 ', 'RDE', datetime(1999,9,10,12,34,56,700000) ),
            (19,  '990910123456.  ', 'RDE', datetime(1999,9,10,12,34,56,000000) ),
            (20,  '990910123456   ', 'RDE', datetime(1999,9,10,12,34,56,000000) ),
            (21,  '9909101234',      'RDE', datetime(1999,9,10,12,34,00,000000) )
        ]

        for (case_no, packed, format_type, answer_time) in time_cases:
            self.subTest(case_no)
            test_time = iod.DateTime_frompacked(packed,format_type)
            self.assertEqual(test_time,answer_time,msg="DateTime_frompacked Case {}  Packed='{}' ({})".format(case_no,packed,format_type))


if __name__ == '__main__':
    main()