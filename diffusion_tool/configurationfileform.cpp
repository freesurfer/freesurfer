#include <QFileInfo>
#include <QDir>
#include <QTextStream>
#include <QFileDialog>
#include <QMessageBox>
#include <QRegularExpression>

#include "configurationfileform.h"

#if (QT_VERSION >= QT_VERSION_CHECK(5, 15, 0))
#define endl Qt::endl
#endif

ConfigurationFileForm::ConfigurationFileForm(QDialog *parent)  : QDialog(parent)
{
    ui.setupUi(this);
    this->setMinimumHeight(859);
    this->setMaximumHeight(859);
    this->setMinimumWidth(783);
    this->setMaximumWidth(783);

    original_dcm = "";
    original_bvec = "";
    original_bval = "";
}

void ConfigurationFileForm::loadDefaults(QString fileName)
{
    QFileInfo file_info(fileName);
    // QTextStream(stdout) << (file_info.dir()).dirName() << " is the name of the directory" << endl;

    QDir parent_directory = file_info.dir();
    QString parent_dir_name = parent_directory.dirName();
    original_directory_path = parent_directory.absolutePath();

    ui.dtrootLineEdit->setText(original_directory_path);
    ui.subjectNameLineEdit->setText(parent_dir_name);
    ui.subjectDcmFileLineEdit->setText(fileName);

    QFileInfoList list = parent_directory.entryInfoList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot);

    ui.subjectBvecFileLineEdit->setText(bvecSeeker(list));
    ui.subjectBvalFileLineEdit->setText(bvalSeeker(list));

    parent_directory.cdUp();
    original_parent_path = parent_directory.absolutePath();

    ui.checkBox->setChecked(true);
    ui.checkBox_2->setChecked(true);
    ui.checkBox_3->setChecked(true);

    ui.doubleSpinBox_2->setValue(0.30);

    ui.checkBox_5->setChecked(true);
    ui.checkBox_7->setChecked(true);

    ui.checkBox_10->setChecked(true);
    ui.checkBox_11->setChecked(true);
    ui.checkBox_12->setChecked(true);
    ui.checkBox_13->setChecked(true);
    ui.checkBox_14->setChecked(true);
    ui.checkBox_15->setChecked(true);
    ui.checkBox_16->setChecked(true);
    ui.checkBox_17->setChecked(true);
    ui.checkBox_18->setChecked(true);
    ui.checkBox_19->setChecked(true);
    ui.checkBox_20->setChecked(true);
    ui.checkBox_21->setChecked(true);
    ui.checkBox_22->setChecked(true);
    ui.checkBox_23->setChecked(true);
    ui.checkBox_24->setChecked(true);
    ui.checkBox_25->setChecked(true);
    ui.checkBox_8->setChecked(true);
    ui.checkBox_9->setChecked(true);

    ui.spinBox->setValue(6);
    ui.spinBox_10->setValue(4);
    ui.spinBox_11->setValue(5);
    ui.spinBox_12->setValue(5);
    ui.spinBox_13->setValue(5);
    ui.spinBox_14->setValue(5);
    ui.spinBox_15->setValue(5);
    ui.spinBox_16->setValue(4);
    ui.spinBox_17->setValue(5);
    ui.spinBox_18->setValue(5);
    ui.spinBox_2->setValue(5);
    ui.spinBox_3->setValue(6);
    ui.spinBox_4->setValue(5);
    ui.spinBox_5->setValue(5);
    ui.spinBox_6->setValue(5);
    ui.spinBox_7->setValue(7);
    ui.spinBox_8->setValue(5);
    ui.spinBox_9->setValue(5);
}

void ConfigurationFileForm::loadDefaults2(QString fileName)
{
    QFileInfo file_info(fileName);
    // QTextStream(stdout) << (file_info.dir()).dirName() << " is the name of the directory" << endl;
    ui.testingLineEdit->setText(fileName);
    QDir temp_directory = file_info.dir();

    QFileInfoList temp_list;
    bool found = false;
    QString end_of_file;

    temp_list = temp_directory.entryInfoList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot);

    // this is a check to see if the ".trk" file exists in the directory
    for (int i = 0; i < temp_list.size(); i++)
    {
        if (found)
        {
            break;
        }

        if (temp_list[i].isFile())
        {
            if (temp_list[i].fileName().length() >= 5)
            {
                end_of_file = temp_list[i].fileName();
                end_of_file.remove(0, end_of_file.length() - 4);

                if (end_of_file == ".trk")
                {
                    found = true;
                    ui.testingLineEdit_2->setText(temp_list[i].absoluteFilePath());
                }

                else
                {
                    QTextStream(stdout) << temp_list[i].absoluteFilePath() << " was just rejected with " << end_of_file << " being the remainder." << endl;
                }
            }
        }
    }

    if (found)
    {
        ui.spinBox_19->setValue(200);
        ui.spinBox_20->setValue(10);
        ui.spinBox_21->setValue(500);

        ui.comboBox->setCurrentIndex(1);
        ui.comboBox_2->setCurrentIndex(2);

        ui.spinBox_22->setValue(1);
    }

    else
    {
        QMessageBox msgBox;
        msgBox.setText("warning: no valid .trk file was found in directory");
        msgBox.exec();
    }
}

void ConfigurationFileForm::on_checkBox_26_stateChanged(int state)
{
    if (state == 2)
    {
        ui.dtrootLineEdit->setText(original_parent_path);

        original_dcm = ui.subjectDcmFileLineEdit->text();
        original_bvec = ui.subjectBvecFileLineEdit->text();
        original_bval = ui.subjectBvalFileLineEdit->text();

        if (!replace_subject_dcm() || !replace_subject_bval() || !replace_subject_bvec())
        {
            QMessageBox msgBox;
            msgBox.setText("warning: subject name was not found in one or more of the dcm, bval, and bvec filepaths");
            msgBox.exec();
        }
    }

    else
    {
        ui.dtrootLineEdit->setText(original_directory_path);

        ui.subjectDcmFileLineEdit->setText(original_dcm);
        ui.subjectBvecFileLineEdit->setText(original_bvec);
        ui.subjectBvalFileLineEdit->setText(original_bval);
    }
}

void ConfigurationFileForm::on_pushButton_clicked(bool checked)
{
    if (!prelimWriteFile())
    {
        return;
    }

    this->setMinimumHeight(1051);
    this->setMaximumHeight(1051);
    QString my_text = "Lines to Run Via Terminal: \n\n";
    my_text += ("1. trac-all -prep -c " + ui.dtrootLineEdit->text() + "/config.tutorial\n\n");
    my_text += ("2. trac-all -bedp -c " + ui.dtrootLineEdit->text() + "/config.tutorial\n\n");
    my_text += ("3. trac-all -path -c " + ui.dtrootLineEdit->text() + "/config.tutorial\n\n");
    ui.textEdit->setText(my_text);
}

void ConfigurationFileForm::on_pushButton_3_clicked(bool checked)
{
    QString full_name = QFileDialog::getOpenFileName(this);
    QFileInfo full_path(full_name);


    if (validDcm(full_path))
    {
        // QTextStream(stdout) << "the file root name is " << file_root << '\n';
        loadDefaults(full_name);
    }

    else
    {
        QMessageBox msgBox;
        msgBox.setText("Must select a DCOM file.");
        msgBox.exec();
    }
}

void ConfigurationFileForm::on_pushButton_2_clicked(bool checked)
{
    QFileInfo dcm_info(ui.testingLineEdit->text());
    QFileInfo trk_info(ui.testingLineEdit_2->text());
    QDir output_dir(ui.outputDirectoryLineEdit->text());

    if (!validDcm(dcm_info))
    {
        QMessageBox msgBox;
        msgBox.setText("Must select a valid DCOM file.");
        msgBox.exec();
        return;
    }

    /*else if (!validTrk(trk_info))
    {
        QMessageBox msgBox;
        msgBox.setText("Must select a valid trk file.");
        msgBox.exec();
        return;
    }*/

    else if (!output_dir.exists())
    {
        QMessageBox msgBox;
        msgBox.setText("Must select a valid output directory.");
        msgBox.exec();
        return;
    }

    else
    {
        this->setMinimumHeight(1051);
        this->setMaximumHeight(1051);
        QString my_text = "Lines to Run Via Terminal: \n\n";
        my_text += ("1. AnatomiCuts -s " + ui.testingLineEdit->text() + " -f " + ui.testingLineEdit_2->text() + " -l " + ui.comboBox_2->currentText()[0] + " -c " + QString::number(ui.spinBox_19->value()) + " -n " + QString::number(ui.spinBox_20->value()) + " -e " + QString::number(ui.spinBox_21->value()) + " -labels -o " + ui.outputDirectoryLineEdit->text() + "\n\n");
        ui.textEdit->setText(my_text);
    }

}

void ConfigurationFileForm::on_pushButton_4_clicked(bool checked)
{
    QString full_name = QFileDialog::getOpenFileName(this);
    QFileInfo full_path(full_name);
    QString file_name = full_path.fileName();

    if (validDcm(full_path))
    {
        loadDefaults2(full_name);
    }

    else
    {
        QMessageBox msgBox;
        msgBox.setText("Must seleerwerwect a valid DCOM file.");
        msgBox.exec();
    }
}

void ConfigurationFileForm::on_buttonBox_2_rejected()
{
    this->close();
}

void ConfigurationFileForm::on_pushButton_5_clicked(bool checked)
{
    this->close();
}

bool ConfigurationFileForm::prelimWriteFile()
{
    // if we're trying to write configuration files for all subjects...
    if (ui.checkBox_26->isChecked())
    {
        QString main_directory = ui.subjectDcmFileLineEdit->text();

        int index = main_directory.indexOf("*");
        main_directory.remove(index, main_directory.size() - index);

        QDir parent_directory(main_directory);

        length_of_root_directory = main_directory.size();

        matcher(parent_directory, ui.subjectDcmFileLineEdit->text(), dcms);
        matcher(parent_directory, ui.subjectBvalFileLineEdit->text(), bvals);
        matcher(parent_directory, ui.subjectBvecFileLineEdit->text(), bvecs);

        writeFile();
        return true;
    }

    // if we're only working with the subject whose file was selected...
    else
    {
        QFileInfo dcm_info(ui.subjectDcmFileLineEdit->text());
        QFileInfo bval_info(ui.subjectBvalFileLineEdit->text());
        QFileInfo bvec_info(ui.subjectBvecFileLineEdit->text());

        //checks if the dcm file exists and if its extension makes sense
        if (!validDcm(dcm_info))
        {
            QMessageBox msgBox;
            msgBox.setText("Must select a valid DCOM file.");
            msgBox.exec();
            return false;
        }

        else if (!validBval(bval_info))
        {
            QMessageBox msgBox;
            msgBox.setText("Must select a valid bval file.");
            msgBox.exec();
            return false;
        }

        else if (!validBvec(bvec_info))
        {
            QMessageBox msgBox;
            msgBox.setText("Must select a valid bvec file.");
            msgBox.exec();
            return false;
        }

        else
        {
            subjects.append(ui.subjectNameLineEdit->text());
            dcms.append(ui.subjectDcmFileLineEdit->text());
            bvals.append(ui.subjectBvalFileLineEdit->text());
            bvecs.append(ui.subjectBvecFileLineEdit->text());
            writeFile();
            return true;
        }
    }
}

QString ConfigurationFileForm::bvecSeeker(QFileInfoList list)
{
    QFileInfo current;

    for (int i = 0; i < list.size(); ++i)
    {
        current = list[i];
        if (current.isFile())
        {
            if (current.fileName() == file_root + ".bvecs" || current.fileName() == file_root + ".voxel_space.bvals")
            {
                return current.filePath();
            }
        }
    }

    return "";
}

QString ConfigurationFileForm::bvalSeeker(QFileInfoList list)
{
    QFileInfo current;

    for (int i = 0; i < list.size(); ++i)
    {
        current = list[i];
        if (current.isFile())
        {
            if (current.fileName() == file_root + ".bvals"  || current.fileName() == file_root + ".voxel_space.bvecs")
            {
                return current.filePath();
            }
        }
    }

    return "";
}

void ConfigurationFileForm::writeFile()
{
    QFile file(ui.dtrootLineEdit->text() + "/config.tutorial");
    // QTextStream(stdout) << "about to start writing the file" << endl;

    if (file.open(QIODevice::ReadWrite))
    {
        QTextStream stream(&file);
        stream << "# FreeSurfer SUBJECTS_DIR" << endl << "# T1 images and FreeSurfer segmentations are expected to be found here" << endl << "#" << endl << "setenv SUBJECTS_DIR " << original_parent_path << endl;
        stream << endl << "# Output directory where trac-all results will be saved" << endl << "# Default: Same as SUBJECTS_DIR" << endl << "#" << endl << "set dtroot = " << ui.dtrootLineEdit->text() << endl;

        // this is the section that sets the subject list
        stream << endl << "# Subject IDs" << endl << "#" << endl << "set subjlist = ( ";
        int count = subjects.length();

        for (int i = 0; i < count; i++)
        {
            if (i == 0)
            {
                stream << subjects[i];
            }

            else
            {
                stream << "\t\t\t\t" << subjects[i];
            }

            if (i == count - 1)
            {
                stream << " )" << endl;
            }

            else
            {
                stream << " \\" << endl;
            }
        }


        // this is the section that provides the .dcm files
        stream << endl << "# Input diffusion DICOMs" << endl << "#" << endl << "set dcmlist = ( ";

        count = dcms.length();

        for (int i = 0; i < count; i++)
        {
            if (i == 0)
            {
                stream << dcms[i];
            }

            else
            {
                stream << "\t\t\t\t" << dcms[i];
            }

              if (i == count - 1)
            {
                stream << " )" << endl;
            }

            else
            {
                stream << " \\" << endl;
            }
        }



        // this is the section that provides the bvec files
        stream << endl << "# Diffusion gradient tables (if there is a different one for each scan)" << endl << "#" << endl << "set bveclist = ( ";

        count = bvecs.length();

        for (int i = 0; i < count; i++)
        {
            if (i == 0)
            {
                stream << bvecs[i];
            }

            else
            {
                stream << "\t\t\t\t" << bvecs[i];
            }

              if (i == count - 1)
            {
                stream << " )" << endl;
            }

            else
            {
                stream << " \\" << endl;
            }
        }



        // this is the section that provides the bval files
        stream << endl << "# Diffusion b-value tables (if there is a different one for each scan)" << endl << "#" << endl << "set bvallist = ( ";

        count = bvals.length();

        for (int i = 0; i < count; i++)
        {
            if (i == 0)
            {
                stream << bvals[i];
            }

            else
            {
                stream << "\t\t\t\t" << bvals[i];
            }

              if (i == count - 1)
            {
                stream << " )" << endl;
            }

            else
            {
                stream << " \\" << endl;
            }
        }



        stream << endl << "set runlist = (";

        for (int i = 0; i < count; i++)
        {
            stream << " " << i + 1;
        }

        stream << " )" << endl;

        stream << endl << "set doeddy = " << ui.checkBox->isChecked() << endl << endl << "set dorotbvecs = " << ui.checkBox_2->isChecked() << endl << endl << "set usemaskanat = " << ui.checkBox_3->isChecked() << endl;
        stream << endl << "set thrbet = " << ui.doubleSpinBox_2->value() << endl << endl << "set doregflt = " << ui.checkBox_4->isChecked() << endl << "set doregbbr = " << ui.checkBox_5->isChecked() << endl << endl << "set doregmni = " << ui.checkBox_7->isChecked() << endl << "set doregcvs = " << ui.checkBox_6->isChecked() << endl;

        QString pathlist = "set pathlist = ( ";

        if (ui.checkBox_9->isChecked())
        {
            pathlist += "lh.cst_AS ";
        }

        if (ui.checkBox_8->isChecked())
        {
            pathlist += "rh.cst_AS ";
        }

        pathlist += "\\ \n \t\t\t\t";

        if (ui.checkBox_13->isChecked())
        {
            pathlist += "lh.ilf_AS ";
        }

        if (ui.checkBox_11->isChecked())
        {
            pathlist += "rh.ilf_AS ";
        }

        pathlist += "\\ \n \t\t\t\t";

        if (ui.checkBox_14->isChecked())
        {
            pathlist += "lh.unc_AS ";
        }

        if (ui.checkBox_15->isChecked())
        {
            pathlist += "rh.unc_AS ";
        }

        pathlist += "\\ \n \t\t\t\t";

        if (ui.checkBox_12->isChecked())
        {
            pathlist += "fmajor_PP ";
        }

        if (ui.checkBox_10->isChecked())
        {
            pathlist += "fminor_PP ";
        }

        pathlist += "\\ \n \t\t\t\t";

        if (ui.checkBox_16->isChecked())
        {
            pathlist += "lh.atr_PP ";
        }

        if (ui.checkBox_18->isChecked())
        {
            pathlist += "rh.atr_PP ";
        }

        pathlist += "\\ \n \t\t\t\t";

        if (ui.checkBox_17->isChecked())
        {
            pathlist += "lh.cab_PP ";
        }

        if (ui.checkBox_21->isChecked())
        {
            pathlist += "rh.cab_PP ";
        }

        pathlist += "\\ \n \t\t\t\t";

        if (ui.checkBox_19->isChecked())
        {
            pathlist += "lh.ccg_PP ";
        }

        if (ui.checkBox_22->isChecked())
        {
            pathlist += "rh.ccg_PP ";
        }

        pathlist += "\\ \n \t\t\t\t";

        if (ui.checkBox_20->isChecked())
        {
            pathlist += "lh.slfp_PP ";
        }

        if (ui.checkBox_24->isChecked())
        {
            pathlist += "rh.slfp_PP ";
        }

        pathlist += "\\ \n \t\t\t\t";

        if (ui.checkBox_23->isChecked())
        {
            pathlist += "lh.slft_PP ";
        }

        if (ui.checkBox_25->isChecked())
        {
            pathlist += "rh.slft_PP ";
        }

        pathlist += ")";
        stream << endl << pathlist << endl;

        QString ncpts = "set ncpts = ( ";

        if (ui.checkBox_9->isChecked())
        {
            ncpts += QString::number(ui.spinBox->value());
            ncpts += " ";
        }

        if (ui.checkBox_8->isChecked())
        {
            ncpts += QString::number(ui.spinBox_3->value());
            ncpts += " ";
        }

        if (ui.checkBox_13->isChecked())
        {
            ncpts += QString::number(ui.spinBox_2->value());
            ncpts += " ";
        }

        if (ui.checkBox_11->isChecked())
        {
            ncpts += QString::number(ui.spinBox_4->value());
            ncpts += " ";
        }

        if (ui.checkBox_14->isChecked())
        {
            ncpts += QString::number(ui.spinBox_5->value());
            ncpts += " ";
        }

        if (ui.checkBox_15->isChecked())
        {
            ncpts += QString::number(ui.spinBox_6->value());
            ncpts += " ";
        }

        if (ui.checkBox_12->isChecked())
        {
            ncpts += QString::number(ui.spinBox_7->value());
            ncpts += " ";
        }

        if (ui.checkBox_10->isChecked())
        {
            ncpts += QString::number(ui.spinBox_13->value());
            ncpts += " ";
        }

        if (ui.checkBox_16->isChecked())
        {
            ncpts += QString::number(ui.spinBox_8->value());
            ncpts += " ";
        }

        if (ui.checkBox_18->isChecked())
        {
            ncpts += QString::number(ui.spinBox_14->value());
            ncpts += " ";
        }

        if (ui.checkBox_17->isChecked())
        {
            ncpts += QString::number(ui.spinBox_9->value());
            ncpts += " ";
        }

        if (ui.checkBox_21->isChecked())
        {
            ncpts += QString::number(ui.spinBox_15->value());
            ncpts += " ";
        }

        if (ui.checkBox_19->isChecked())
        {
            ncpts += QString::number(ui.spinBox_10->value());
            ncpts += " ";
        }

        if (ui.checkBox_22->isChecked())
        {
            ncpts += QString::number(ui.spinBox_16->value());
            ncpts += " ";
        }

        if (ui.checkBox_20->isChecked())
        {
            ncpts += QString::number(ui.spinBox_11->value());
            ncpts += " ";
        }

        if (ui.checkBox_24->isChecked())
        {
            ncpts += QString::number(ui.spinBox_17->value());
            ncpts += " ";
        }

        if (ui.checkBox_23->isChecked())
        {
            ncpts += QString::number(ui.spinBox_12->value());
            ncpts += " ";
        }

        if (ui.checkBox_25->isChecked())
        {
            ncpts += QString::number(ui.spinBox_18->value());
            ncpts += " ";
        }

        ncpts += ")";

        stream << endl << ncpts << endl;
    }

    else
    {
        QTextStream(stdout) << "could not write configuration file" << endl;
    }
}

bool ConfigurationFileForm::validDcm(QFileInfo file)
{
    if (!file.exists())
    {
        return false;
    }

    QString name = file.fileName();
    QString end_of_file;
    int length = name.length();

    if (length <= 4)
    {
        return false;
    }

    if (length > 7)
    {
        end_of_file = name;
        end_of_file.remove(0, length - 7);

        // checking to see if ".nii.gz" is extension
        if (end_of_file == ".nii.gz")
        {
            name.remove(length - 7, length);
            file_root = name;
            return true;
        }
    }

    end_of_file = name;
    end_of_file.remove(0, length - 4);

    // checking to see if ".dcm" is extension
    if (end_of_file == ".dcm")
    {
        name.remove(length - 4, length);
        file_root = name;
        return true;
    }

    else
    {
        return false;
    }
}

bool ConfigurationFileForm::replace_subject_dcm()
{
    QString full_line = ui.subjectDcmFileLineEdit->text();
    QString query_string = ui.subjectNameLineEdit->text();
    int query_length = query_string.size();

    int pos = 0;
    int index = full_line.indexOf(query_string, pos);

    if (index == -1)
    {
        return false;
    }

    while (index != -1)
    {
        full_line.replace(index, query_length, "*");
        pos += (index = 1);
        index = full_line.indexOf(query_string, pos);
    }

    ui.subjectDcmFileLineEdit->setText(full_line);

    return true;
}

bool ConfigurationFileForm::replace_subject_bvec()
{
    QString full_line = ui.subjectBvecFileLineEdit->text();
    QString query_string = ui.subjectNameLineEdit->text();
    int query_length = query_string.size();

    int pos = 0;
    int index = full_line.indexOf(query_string, pos);

    if (index == -1)
    {
        return false;
    }

    while (index != -1)
    {
        full_line.replace(index, query_length, "*");
        pos += (index = 1);
        index = full_line.indexOf(query_string, pos);
    }

    ui.subjectBvecFileLineEdit->setText(full_line);

    return true;
}

bool ConfigurationFileForm::replace_subject_bval()
{
    QString full_line = ui.subjectBvalFileLineEdit->text();
    QString query_string = ui.subjectNameLineEdit->text();
    int query_length = query_string.size();

    int pos = 0;
    int index = full_line.indexOf(query_string, pos);

    if (index == -1)
    {
        return false;
    }

    while (index != -1)
    {
        full_line.replace(index, query_length, "*");
        pos += (index = 1);
        index = full_line.indexOf(query_string, pos);
    }

    ui.subjectBvalFileLineEdit->setText(full_line);
    return true;
}

void ConfigurationFileForm::matcher(QDir directory, QString pattern, QVector<QString>& vector)
{
    QFileInfoList directory_contents = directory.entryInfoList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot);
    // QTextStream(stdout) << "about to peep " << directory.absolutePath() << " and its length is " << directory_contents.length() << endl;
    QDir* temp = NULL;
    QString wc_exp = QRegularExpression::wildcardToRegularExpression(pattern);
    QRegularExpression rx(QRegularExpression::anchoredPattern(wc_exp),
                          QRegularExpression::CaseInsensitiveOption);
    QString subject_name;
    int temp_index;

    for (int i = 0; i < directory_contents.length(); i++)
    {
        if (directory_contents[i].isDir())
        {
            temp = new QDir(directory_contents[i].absoluteFilePath());
            // QTextStream(stdout) << "the name of this directory is " << temp->absolutePath() << endl;
            matcher(*temp, pattern, vector);
            delete temp;
            temp = NULL;
        }

        else
        {
            // QTextStream(stdout) << "about to peep the contents of " << directory_contents[i].absoluteFilePath() << endl;
            if (rx.match(directory_contents[i].absoluteFilePath()).hasMatch())
            {
                vector.append(directory_contents[i].absoluteFilePath());

                if (vector == dcms)
                {
                    // removing the stuff that comes after the subject name
                    subject_name = directory_contents[i].absoluteFilePath();
                    temp_index = subject_name.indexOf("/", length_of_root_directory);
                    subject_name.remove(temp_index, subject_name.size() - temp_index);

                    // removing the stuff that comes before the subject name
                    temp_index = subject_name.lastIndexOf("/");
                    subject_name.remove(0, temp_index + 1);

                    // adding the subject name to the right vector
                    subjects.append(subject_name);
                }
            }
        }
    }
}

bool validBval(QFileInfo file)
{
    if (!file.exists())
    {
        return false;
    }

    QString name = file.fileName();
    QString end_of_file;
    int length = name.length();

    if (length <= 6)
    {
        return false;
    }

    if (length > 18)
    {
        end_of_file = name;
        end_of_file.remove(0, length - 18);

        // checking to see if ".voxel_space.bvals" is extension
        if (end_of_file == ".voxel_space.bvals")
        {
            return true;
        }
    }

    end_of_file = name;
    end_of_file.remove(0, length - 6);

    // checking to see if ".bvals"
    if (end_of_file == ".bvals")
    {
        return true;
    }

    else
    {
        return false;
    }
}

bool validBvec(QFileInfo file)
{
    if (!file.exists())
    {
        return false;
    }

    QString name = file.fileName();
    QString end_of_file;
    int length = name.length();

    if (length <= 6)
    {
        return false;
    }

    if (length > 18)
    {
        end_of_file = name;
        end_of_file.remove(0, length - 18);

        // checking to see if ".voxel_space.bvecs" is extension
        if (end_of_file == ".voxel_space.bvecs")
        {
            return true;
        }
    }

    end_of_file = name;
    end_of_file.remove(0, length - 6);

    // checking to see if ".bvecs"
    if (end_of_file == ".bvecs")
    {
        return true;
    }

    else
    {
        return false;
    }
}

bool validTrk(QFileInfo file)
{
    if (!file.exists())
    {
        return false;
    }

    QString name = file.fileName();
    QString end_of_file;
    int length = name.length();

    if (length <= 4)
    {
        return false;
    }

    end_of_file = name;
    end_of_file.remove(0, length - 4);

    // checking to see if ".trk"
    if (end_of_file == "trk")
    {
        return true;
    }

    else
    {
        return false;
    }
}
