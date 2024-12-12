from .StandardAmpliconWorkflow import StandardAmpliconWorkflow
from .StandardMetagenomicWorkflow import StandardMetagenomicWorkflow
from .StandardMetatranscriptomicWorkflow import \
    StandardMetatranscriptomicWorkflow
from .TellseqMetagenomicWorkflow import TellSeqMetagenomicWorkflow
from sequence_processing_pipeline.Pipeline import Pipeline
from metapool import load_sample_sheet
from .Assays import METAOMIC_ASSAY_NAMES, ASSAY_NAME_AMPLICON
from .Protocol import PROTOCOL_NAME_ILLUMINA, PROTOCOL_NAME_TELLSEQ
from .Workflows import WorkflowError


class WorkflowFactory():
    WORKFLOWS = [StandardMetagenomicWorkflow,
                 StandardMetatranscriptomicWorkflow,
                 StandardAmpliconWorkflow,
                 TellSeqMetagenomicWorkflow]

    ST_TO_IN_MAP = {PROTOCOL_NAME_ILLUMINA: ['standard_metag',
                                             'standard_metat',
                                             'absquant_metag',
                                             'absquant_metat'],
                    PROTOCOL_NAME_TELLSEQ: ['tellseq_metag',
                                            'tellseq_absquant']}

    @classmethod
    def _get_instrument_type(cls, sheet):
        for instrument_type in cls.ST_TO_IN_MAP:
            if sheet.Header['SheetType'] in cls.ST_TO_IN_MAP[instrument_type]:
                return instrument_type

    @classmethod
    def generate_workflow(cls, **kwargs):
        msg = "kwargs must not be None and must define 'uif_path'"

        if not kwargs:
            # if kwargs is None or {}, raise an Error
            raise ValueError(msg)

        if 'uif_path' not in kwargs:
            raise ValueError(msg)

        if Pipeline.is_sample_sheet(kwargs['uif_path']):
            # NB: The Pipeline() determines an input-file is a sample-sheet
            # if the first line begins with "[Header]" followed by any number
            # of ','. A file that begins this way but fails to load
            # successfully because of an undefined SheetType and/or
            # SheetVersion will raise a ValueError() here, w/the message
            # "'{sheet}' doesn't appear to be a valid sample-sheet."

            sheet = load_sample_sheet(kwargs['uif_path'])

            # if we do not validate the sample-sheet now, it will be validated
            # downstream when we attempt to instantiate a Workflow(), which in
            # turn will attempt to instantiate a Pipeline(), which will load
            # and validate the sample-sheet on its own. This is an early
            # abort. Expect the user/caller to diagnose the sample-sheet in a
            # notebook or by other means.
            if sheet.validate_and_scrub_sample_sheet():
                assay_type = sheet.Header['Assay']
                if assay_type not in METAOMIC_ASSAY_NAMES:
                    # NB: This Error is not likely to be raised unless an
                    # assay type is defined in metapool but not in Assays.
                    raise WorkflowError("Can't determine workflow from assay "
                                        "type: %s" % assay_type)
                instrument_type = cls._get_instrument_type(sheet)
            else:
                raise WorkflowError(f"'{kwargs['uif_path']} doesn't appear to "
                                    "be a valid sample-sheet.")
        elif Pipeline.is_mapping_file(kwargs['uif_path']):
            # if file is readable as a basic TSV and contains all the required
            # headers, then treat this as a mapping file, even if it's an
            # invalid one.
            assay_type = ASSAY_NAME_AMPLICON
            # for Amplicon runs, the lane_number is always one, even if the
            # user supplies another value in the UI.
            kwargs['lane_number'] = 1
            # NB: For now, let's assume all Amplicon runs are Illumina, since
            # the entire Amplicon pipeline assumes as much.
            instrument_type = 'Illumina'
        else:
            raise ValueError("Your uploaded file doesn't appear to be a "
                             "sample-sheet or a mapping-file.")

        for workflow in WorkflowFactory.WORKFLOWS:
            if workflow.assay_type == assay_type:
                if workflow.protocol_type == instrument_type:
                    # return instantiated workflow object
                    return workflow(**kwargs)

        # This Error will only be raised if a sample-sheet passes metapool's
        # validation method but a Workflow() for its instrument-type and
        # assay-type doesn't exist.
        raise ValueError(f"Assay type '{assay_type}' and Instrument type "
                         f"'{instrument_type}' did not match any known "
                         "workflow configuration")
