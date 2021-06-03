SELECT DISTINCT 
                    CASE IDMapping.is_canonical
                        WHEN 1
                        THEN IDMapping.Entry
                        ELSE IDMapping.isoform
                    END UniProt,Mutation.Ref, Mutation.Pos, Mutation.Alt,
                    Mutation.Pos - ResidueMappingRange.unp_beg + ResidueMappingRange.pdb_beg AS residue_number,
                    Mutation.Pos - ResidueMappingRange.unp_beg + ResidueMappingRange.auth_pdb_beg AS author_residue_number,
                    ResidueMappingRange.author_insertion_code,
                    ResidueMappingRange.observed_ratio,
                    ResidueMappingRange.pdb_id,
                    ResidueMappingRange.entity_id,
                    ResidueMappingRange.struct_asym_id,
                    ResidueMappingRange.chain_id,
                    (CASE
                        WHEN ResidueMappingRange.residue_name == ''
                        THEN (SELECT three_letter_code FROM AAThree2one WHERE one_letter_code == UniProtSeq.Ref)
                        ELSE ResidueMappingRange.residue_name
                        END
                    ) AS residue_name,
                    UniProtSeq.Ref AS unp_one_letter_code,
                    (CASE
                        WHEN ResidueMappingRange.conflict_code IS NULL
                        THEN UniProtSeq.Ref
                        ELSE (
                            CASE
                            WHEN ResidueMappingRange.residue_name IN (SELECT three_letter_code FROM AAThree2one)
                            THEN (SELECT one_letter_code FROM AAThree2one WHERE three_letter_code == ResidueMappingRange.residue_name)
                            ELSE 'X'
                            END
                        )
                        END
                    ) AS pdb_one_letter_code
        FROM Mutation,ResidueMappingRange
            INNER JOIN IDMapping ON Mutation.ftId = IDMapping.ftId
            INNER JOIN UniProtSeq ON UniProtSeq.isoform = IDMapping.isoform 
                                AND UniProtSeq.Pos = Mutation.Pos 
        WHERE ResidueMappingRange.UniProt = UniProt
        AND Mutation.Pos >= ResidueMappingRange.unp_beg
        AND Mutation.Pos <= ResidueMappingRange.unp_end
        ;